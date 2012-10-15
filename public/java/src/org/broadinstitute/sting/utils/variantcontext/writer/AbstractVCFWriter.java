package org.broadinstitute.sting.utils.variantcontext.writer;

import java.io.File;
import java.io.OutputStream;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import net.sf.samtools.SAMSequenceDictionary;

public abstract class AbstractVCFWriter
	extends IndexingVariantContextWriter
	{
    // the VCF header we're storing
	protected VCFHeader mHeader = null;

	protected IntGenotypeFieldAccessors intGenotypeFieldAccessors = new IntGenotypeFieldAccessors();

    // should we write genotypes or just sites?
    final protected boolean doNotWriteGenotypes;
    final protected boolean allowMissingFieldsInHeader;

	
    protected AbstractVCFWriter(
    		final File location,
    		final OutputStream output,
    		final SAMSequenceDictionary refDict,
            final boolean enableOnTheFlyIndexing,
            boolean doNotWriteGenotypes,
            final boolean allowMissingFieldsInHeader
            )
    	{
    	super(writerName(location, output), location, output, refDict, enableOnTheFlyIndexing);
        this.doNotWriteGenotypes = doNotWriteGenotypes;
        this.allowMissingFieldsInHeader = allowMissingFieldsInHeader;

    	}
    
    protected VCFHeader getVCFHeader()
    	{
    	return this.mHeader;
    	}
    
   protected static Map<Allele, String> buildAlleleMap(final VariantContext vc) {
        final Map<Allele, String> alleleMap = new HashMap<Allele, String>(vc.getAlleles().size()+1);
        alleleMap.put(Allele.NO_CALL, VCFConstants.EMPTY_ALLELE); // convenience for lookup

        final List<Allele> alleles = vc.getAlleles();
        for ( int i = 0; i < alleles.size(); i++ ) {
            alleleMap.put(alleles.get(i), String.valueOf(i));
        }

        return alleleMap;
    }

   
   private static final String QUAL_FORMAT_STRING = "%.2f";
   private static final String QUAL_FORMAT_EXTENSION_TO_TRIM = ".00";

   protected String formatQualValue(double qual) {
       String s = String.format(QUAL_FORMAT_STRING, qual);
       if ( s.endsWith(QUAL_FORMAT_EXTENSION_TO_TRIM) )
           s = s.substring(0, s.length() - QUAL_FORMAT_EXTENSION_TO_TRIM.length());
       return s;
   }

   public static final void missingSampleError(final VariantContext vc, final VCFHeader header) {
       final List<String> badSampleNames = new ArrayList<String>();
       for ( final String x : header.getGenotypeSamples() )
           if ( ! vc.hasGenotype(x) ) badSampleNames.add(x);
       throw new ReviewedStingException("BUG: we now require all samples in VCFheader to have genotype objects.  Missing samples are " + Utils.join(",", badSampleNames));
   }
   
   protected boolean isMissingValue(String s) {
       // we need to deal with the case that it's a list of missing values
       return (countOccurrences(VCFConstants.MISSING_VALUE_v4.charAt(0), s) + countOccurrences(',', s) == s.length());
   }

   
   /**
    * Takes a double value and pretty prints it to a String for display
    *
    * Large doubles => gets %.2f style formatting
    * Doubles < 1 / 10 but > 1/100 </>=> get %.3f style formatting
    * Double < 1/100 => %.3e formatting
    * @param d
    * @return
    */
   public static final String formatVCFDouble(final double d) {
       String format;
       if ( d < 1 ) {
           if ( d < 0.01 ) {
               if ( Math.abs(d) >= 1e-20 )
                   format = "%.3e";
               else {
                   // return a zero format
                   return "0.00";
               }
           } else {
               format = "%.3f";
           }
       } else {
           format = "%.2f";
       }

       return String.format(format, d);
   }

   public static String formatVCFField(Object val) {
       String result;
       if ( val == null )
           result = VCFConstants.MISSING_VALUE_v4;
       else if ( val instanceof Double )
           result = formatVCFDouble((Double) val);
       else if ( val instanceof Boolean )
           result = (Boolean)val ? "" : null; // empty string for true, null for false
       else if ( val instanceof List ) {
           result = formatVCFField(((List)val).toArray());
       } else if ( val.getClass().isArray() ) {
           int length = Array.getLength(val);
           if ( length == 0 )
               return formatVCFField(null);
           StringBuffer sb = new StringBuffer(formatVCFField(Array.get(val, 0)));
           for ( int i = 1; i < length; i++) {
               sb.append(",");
               sb.append(formatVCFField(Array.get(val, i)));
           }
           result = sb.toString();
       } else
           result = val.toString();

       return result;
   }

   /**
    * Determine which genotype fields are in use in the genotypes in VC
    * @param vc
    * @return an ordered list of genotype fields in use in VC.  If vc has genotypes this will always include GT first
    */
   public static List<String> calcVCFGenotypeKeys(final VariantContext vc, final VCFHeader header) {
       Set<String> keys = new HashSet<String>();

       boolean sawGoodGT = false;
       boolean sawGoodQual = false;
       boolean sawGenotypeFilter = false;
       boolean sawDP = false;
       boolean sawAD = false;
       boolean sawPL = false;
       for ( final Genotype g : vc.getGenotypes() ) {
           keys.addAll(g.getExtendedAttributes().keySet());
           if ( g.isAvailable() ) sawGoodGT = true;
           if ( g.hasGQ() ) sawGoodQual = true;
           if ( g.hasDP() ) sawDP = true;
           if ( g.hasAD() ) sawAD = true;
           if ( g.hasPL() ) sawPL = true;
           if (g.isFiltered()) sawGenotypeFilter = true;
       }

       if ( sawGoodQual ) keys.add(VCFConstants.GENOTYPE_QUALITY_KEY);
       if ( sawDP ) keys.add(VCFConstants.DEPTH_KEY);
       if ( sawAD ) keys.add(VCFConstants.GENOTYPE_ALLELE_DEPTHS);
       if ( sawPL ) keys.add(VCFConstants.GENOTYPE_PL_KEY);
       if ( sawGenotypeFilter ) keys.add(VCFConstants.GENOTYPE_FILTER_KEY);

       List<String> sortedList = ParsingUtils.sortList(new ArrayList<String>(keys));

       // make sure the GT is first
       if ( sawGoodGT ) {
           List<String> newList = new ArrayList<String>(sortedList.size()+1);
           newList.add(VCFConstants.GENOTYPE_KEY);
           newList.addAll(sortedList);
           sortedList = newList;
       }

       if ( sortedList.isEmpty() && header.hasGenotypingData() ) {
           // this needs to be done in case all samples are no-calls
           return Collections.singletonList(VCFConstants.GENOTYPE_KEY);
       } else {
           return sortedList;
       }
   }


   private static int countOccurrences(char c, String s) {
          int count = 0;
          for (int i = 0; i < s.length(); i++) {
              count += s.charAt(i) == c ? 1 : 0;
          }
          return count;
   }


	}
