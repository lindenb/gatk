/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.variantcontext.writer;

import net.sf.samtools.SAMSequenceDictionary;
import org.broad.tribble.TribbleException;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.io.*;
import java.util.*;

/**
 * this class writes VCF files
 */
class VCFWriter extends AbstractVCFWriter {
    private final static String VERSION_LINE = VCFHeader.METADATA_INDICATOR + VCFHeaderVersion.VCF4_1.getFormatString() + "=" + VCFHeaderVersion.VCF4_1.getVersionString();

    // the print stream we're writing to
    final protected BufferedWriter mWriter;




    public VCFWriter(final File location, final OutputStream output, final SAMSequenceDictionary refDict,
                     final boolean enableOnTheFlyIndexing,
                     boolean doNotWriteGenotypes,
                     final boolean allowMissingFieldsInHeader )
    	{
        super(location, output, refDict, enableOnTheFlyIndexing,doNotWriteGenotypes,allowMissingFieldsInHeader);
        mWriter = new BufferedWriter(new OutputStreamWriter(getOutputStream())); // todo -- fix buffer size
    	}

    // --------------------------------------------------------------------------------
    //
    // VCFWriter interface functions
    //
    // --------------------------------------------------------------------------------

    @Override
    public void writeHeader(VCFHeader header) {
        // note we need to update the mHeader object after this call because they header
        // may have genotypes trimmed out of it, if doNotWriteGenotypes is true
        super.mHeader = writeHeader(header, mWriter, doNotWriteGenotypes, getVersionLine(), getStreamName());
    	}

    public static final String getVersionLine() {
        return VERSION_LINE;
    }

    public static VCFHeader writeHeader(VCFHeader header,
                                        final Writer writer,
                                        final boolean doNotWriteGenotypes,
                                        final String versionLine,
                                        final String streamNameForError) {
        header = doNotWriteGenotypes ? new VCFHeader(header.getMetaDataInSortedOrder()) : header;
        
        try {
            // the file format field needs to be written first
            writer.write(versionLine + "\n");

            for ( VCFHeaderLine line : header.getMetaDataInSortedOrder() ) {
                if ( VCFHeaderVersion.isFormatString(line.getKey()) )
                    continue;

                writer.write(VCFHeader.METADATA_INDICATOR);
                writer.write(line.toString());
                writer.write("\n");
            }

            // write out the column line
            writer.write(VCFHeader.HEADER_INDICATOR);
            boolean isFirst = true;
            for ( VCFHeader.HEADER_FIELDS field : header.getHeaderFields() ) {
                if ( isFirst )
                    isFirst = false; // don't write out a field separator
                else
                    writer.write(VCFConstants.FIELD_SEPARATOR);
                writer.write(field.toString());
            }

            if ( header.hasGenotypingData() ) {
                writer.write(VCFConstants.FIELD_SEPARATOR);
                writer.write("FORMAT");
                for ( String sample : header.getGenotypeSamples() ) {
                    writer.write(VCFConstants.FIELD_SEPARATOR);
                    writer.write(sample);
                }
            }

            writer.write("\n");
            writer.flush();  // necessary so that writing to an output stream will work
        }
        catch (IOException e) {
            throw new ReviewedStingException("IOException writing the VCF header to " + streamNameForError, e);
        }

        return header;
    }

    /**
     * attempt to close the VCF file
     */
    @Override
    public void close() {
        // try to close the vcf stream
        try {
            mWriter.flush();
            mWriter.close();
        } catch (IOException e) {
            throw new ReviewedStingException("Unable to close " + getStreamName(), e);
        }

        super.close();
    }

    /**
     * add a record to the file
     *
     * @param vc      the Variant Context object
     */
    @Override
    public void add(VariantContext vc) {
        if ( mHeader == null )
            throw new IllegalStateException("The VCF Header must be written before records can be added: " + getStreamName());

        if ( doNotWriteGenotypes )
            vc = new VariantContextBuilder(vc).noGenotypes().make();

        try {
            super.add(vc);

            Map<Allele, String> alleleMap = buildAlleleMap(vc);

            // CHROM
            mWriter.write(vc.getChr());
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // POS
            mWriter.write(String.valueOf(vc.getStart()));
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // ID
            String ID = vc.getID();
            mWriter.write(ID);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // REF
            String refString = vc.getReference().getDisplayString();
            mWriter.write(refString);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // ALT
            if ( vc.isVariant() ) {
                Allele altAllele = vc.getAlternateAllele(0);
                String alt = altAllele.getDisplayString();
                mWriter.write(alt);

                for (int i = 1; i < vc.getAlternateAlleles().size(); i++) {
                    altAllele = vc.getAlternateAllele(i);
                    alt = altAllele.getDisplayString();
                    mWriter.write(",");
                    mWriter.write(alt);
                }
            } else {
                mWriter.write(VCFConstants.EMPTY_ALTERNATE_ALLELE_FIELD);
            }
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // QUAL
            if ( !vc.hasLog10PError() )
                mWriter.write(VCFConstants.MISSING_VALUE_v4);
            else
                mWriter.write(formatQualValue(vc.getPhredScaledQual()));
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // FILTER
            String filters = getFilterString(vc);
            mWriter.write(filters);
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            // INFO
            Map<String, String> infoFields = new TreeMap<String, String>();
            for ( Map.Entry<String, Object> field : vc.getAttributes().entrySet() ) {
                String key = field.getKey();

                if ( ! mHeader.hasInfoLine(key) )
                    fieldIsMissingFromHeaderError(vc, key, "INFO");

                String outputValue = formatVCFField(field.getValue());
                if ( outputValue != null )
                    infoFields.put(key, outputValue);
            }
            writeInfoString(infoFields);

            // FORMAT
            final GenotypesContext gc = vc.getGenotypes();
            if ( gc.isLazyWithData() && ((LazyGenotypesContext)gc).getUnparsedGenotypeData() instanceof String ) {
                mWriter.write(VCFConstants.FIELD_SEPARATOR);
                mWriter.write(((LazyGenotypesContext)gc).getUnparsedGenotypeData().toString());
            } else {
                List<String> genotypeAttributeKeys = calcVCFGenotypeKeys(vc, mHeader);
                if ( ! genotypeAttributeKeys.isEmpty() ) {
                    for ( final String format : genotypeAttributeKeys )
                        if ( ! mHeader.hasFormatLine(format) )
                            fieldIsMissingFromHeaderError(vc, format, "FORMAT");

                    final String genotypeFormatString = ParsingUtils.join(VCFConstants.GENOTYPE_FIELD_SEPARATOR, genotypeAttributeKeys);

                    mWriter.write(VCFConstants.FIELD_SEPARATOR);
                    mWriter.write(genotypeFormatString);

                    addGenotypeData(vc, alleleMap, genotypeAttributeKeys);
                }
            }
            
            mWriter.write("\n");
            mWriter.flush();  // necessary so that writing to an output stream will work
        } catch (IOException e) {
            throw new RuntimeException("Unable to write the VCF object to " + getStreamName());
        }
    }


    // --------------------------------------------------------------------------------
    //
    // implementation functions
    //
    // --------------------------------------------------------------------------------

    private final String getFilterString(final VariantContext vc) {
        if ( vc.isFiltered() ) {
            for ( final String filter : vc.getFilters() )
                if ( ! mHeader.hasFilterLine(filter) )
                    fieldIsMissingFromHeaderError(vc, filter, "FILTER");

            return ParsingUtils.join(";", ParsingUtils.sortList(vc.getFilters()));
        }
        else if ( vc.filtersWereApplied() )
            return VCFConstants.PASSES_FILTERS_v4;
        else
            return VCFConstants.UNFILTERED;
    }


    /**
     * create the info string; assumes that no values are null
     *
     * @param infoFields a map of info fields
     * @throws IOException for writer
     */
    private void writeInfoString(Map<String, String> infoFields) throws IOException {
        if ( infoFields.isEmpty() ) {
            mWriter.write(VCFConstants.EMPTY_INFO_FIELD);
            return;
        }

        boolean isFirst = true;
        for ( Map.Entry<String, String> entry : infoFields.entrySet() ) {
            if ( isFirst )
                isFirst = false;
            else
                mWriter.write(VCFConstants.INFO_FIELD_SEPARATOR);

            String key = entry.getKey();
            mWriter.write(key);

            if ( !entry.getValue().equals("") ) {
                VCFInfoHeaderLine metaData = mHeader.getInfoHeaderLine(key);
                if ( metaData == null || metaData.getCountType() != VCFHeaderLineCount.INTEGER || metaData.getCount() != 0 ) {
                    mWriter.write("=");
                    mWriter.write(entry.getValue());
                }
            }
        }
    }

    /**
     * add the genotype data
     *
     * @param vc                     the variant context
     * @param genotypeFormatKeys  Genotype formatting string
     * @param alleleMap              alleles for this context
     * @throws IOException for writer
     */
    private void addGenotypeData(VariantContext vc, Map<Allele, String> alleleMap, List<String> genotypeFormatKeys)
    throws IOException {
        for ( String sample : mHeader.getGenotypeSamples() ) {
            mWriter.write(VCFConstants.FIELD_SEPARATOR);

            Genotype g = vc.getGenotype(sample);
            if ( g == null ) {
                missingSampleError(vc, mHeader);
            }

            final List<String> attrs = new ArrayList<String>(genotypeFormatKeys.size());
            for ( String field : genotypeFormatKeys ) {
                if ( field.equals(VCFConstants.GENOTYPE_KEY) ) {
                    if ( !g.isAvailable() ) {
                        throw new ReviewedStingException("GTs cannot be missing for some samples if they are available for others in the record");
                    }

                    writeAllele(g.getAllele(0), alleleMap);
                    for (int i = 1; i < g.getPloidy(); i++) {
                        mWriter.write(g.isPhased() ? VCFConstants.PHASED : VCFConstants.UNPHASED);
                        writeAllele(g.getAllele(i), alleleMap);
                    }

                    continue;
                } else {
                    String outputValue;
                    if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY ) ) {
                        outputValue = g.isFiltered() ? g.getFilters() : VCFConstants.PASSES_FILTERS_v4;
                    } else {
                        final IntGenotypeFieldAccessors.Accessor accessor = intGenotypeFieldAccessors.getAccessor(field);
                        if ( accessor != null ) {
                            final int[] intValues = accessor.getValues(g);
                            if ( intValues == null )
                                outputValue = VCFConstants.MISSING_VALUE_v4;
                            else if ( intValues.length == 1 ) // fast path
                                outputValue = Integer.toString(intValues[0]);
                            else {
                                StringBuilder sb = new StringBuilder();
                                sb.append(intValues[0]);
                                for ( int i = 1; i < intValues.length; i++) {
                                    sb.append(",");
                                    sb.append(intValues[i]);
                                }
                                outputValue = sb.toString();
                            }
                        } else {
                            Object val = g.hasExtendedAttribute(field) ? g.getExtendedAttribute(field) : VCFConstants.MISSING_VALUE_v4;

                            VCFFormatHeaderLine metaData = mHeader.getFormatHeaderLine(field);
                            if ( metaData != null ) {
                                int numInFormatField = metaData.getCount(vc);
                                if ( numInFormatField > 1 && val.equals(VCFConstants.MISSING_VALUE_v4) ) {
                                    // If we have a missing field but multiple values are expected, we need to construct a new string with all fields.
                                    // For example, if Number=2, the string has to be ".,."
                                    StringBuilder sb = new StringBuilder(VCFConstants.MISSING_VALUE_v4);
                                    for ( int i = 1; i < numInFormatField; i++ ) {
                                        sb.append(",");
                                        sb.append(VCFConstants.MISSING_VALUE_v4);
                                    }
                                    val = sb.toString();
                                }
                            }

                            // assume that if key is absent, then the given string encoding suffices
                            outputValue = formatVCFField(val);
                        }
                    }

                    if ( outputValue != null )
                        attrs.add(outputValue);
                }
            }

            // strip off trailing missing values
            for (int i = attrs.size()-1; i >= 0; i--) {
                if ( isMissingValue(attrs.get(i)) )
                    attrs.remove(i);
                else
                    break;
            }

            for (int i = 0; i < attrs.size(); i++) {
                if ( i > 0 || genotypeFormatKeys.contains(VCFConstants.GENOTYPE_KEY) )
                    mWriter.write(VCFConstants.GENOTYPE_FIELD_SEPARATOR);
                mWriter.write(attrs.get(i));
            }
        }
    }



    private void writeAllele(Allele allele, Map<Allele, String> alleleMap) throws IOException {
        String encoding = alleleMap.get(allele);
        if ( encoding == null )
            throw new TribbleException.InternalCodecException("Allele " + allele + " is not an allele in the variant context");
        mWriter.write(encoding);
    }

    private final void fieldIsMissingFromHeaderError(final VariantContext vc, final String id, final String field) {
        if ( !allowMissingFieldsInHeader)
            throw new UserException.MalformedVCFHeader("Key " + id + " found in VariantContext field " + field
                    + " at " + vc.getChr() + ":" + vc.getStart()
                    + " but this key isn't defined in the VCFHeader.  The GATK now requires all VCFs to have"
                    + " complete VCF headers by default.  This error can be disabled with the engine argument"
                    + " -U LENIENT_VCF_PROCESSING");
    }
}
