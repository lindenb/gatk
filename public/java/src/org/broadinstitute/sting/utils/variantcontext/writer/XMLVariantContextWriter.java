package org.broadinstitute.sting.utils.variantcontext.writer;

import java.io.File;
import java.io.OutputStream;
import java.util.Map;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.broadinstitute.sting.utils.codecs.vcf.VCFContigHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFilterHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderVersion;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import net.sf.samtools.SAMSequenceDictionary;

public class XMLVariantContextWriter
	extends AbstractVCFWriter
	{
	public final String NS="http://xml.1000genomes.org/";
    // the print stream we're writing to
    final protected XMLStreamWriter mWriter;

	public XMLVariantContextWriter(final File location, final OutputStream output, final SAMSequenceDictionary refDict,
            final boolean enableOnTheFlyIndexing,
            boolean doNotWriteGenotypes,
            final boolean allowMissingFieldsInHeader )
			throws XMLStreamException
		{
		super(location, output, refDict, enableOnTheFlyIndexing,doNotWriteGenotypes,allowMissingFieldsInHeader);
		XMLOutputFactory factory=XMLOutputFactory.newInstance();
		this.mWriter=factory.createXMLStreamWriter(super.getOutputStream());
		}

	protected void writeMetaData(String key,String value)
		throws XMLStreamException
		{
		if(value!=null)
			{
	    	this.mWriter.writeStartElement("metadata");
	    	this.mWriter.writeAttribute("key",key);
	    	this.mWriter.writeCharacters(value);
	    	this.mWriter.writeEndElement();
			}
		else
			{
	    	this.mWriter.writeEmptyElement("metadata");
	    	this.mWriter.writeAttribute("key",key);
			}
		}
	
	@Override
	public void writeHeader(VCFHeader header)
		{
   
    	//
    	
        header = doNotWriteGenotypes ? new VCFHeader(header.getMetaDataInSortedOrder()) : header;
        
        try {
        	this.mWriter.writeStartElement("vcf");
        	this.mWriter.writeAttribute("xmlns", NS);
        	this.mWriter.writeStartElement("head");
        	
        	writeMetaData(
        			VCFHeaderVersion.VCF4_1.getFormatString(),
        			VCFHeaderVersion.VCF4_1.getVersionString()
        			);
        	//INFO
        	this.mWriter.writeStartElement("info-list");
        	for ( VCFInfoHeaderLine line : header.getInfoHeaderLines() )
        	  	{
            	this.mWriter.writeStartElement("info");
            	this.mWriter.writeAttribute("ID",line.getID());
            	this.mWriter.writeAttribute("type",line.getType().name());
            	if(line.isFixedCount()) this.mWriter.writeAttribute("count",String.valueOf(line.getCount()));
            	this.mWriter.writeCharacters(line.getDescription());
            	this.mWriter.writeEndElement();
        	  	}
        	this.mWriter.writeEndElement();
        	
        	//FORMAT
        	this.mWriter.writeStartElement("format-list");
        	for ( VCFFormatHeaderLine line : header.getFormatHeaderLines() )
        	  	{
            	this.mWriter.writeStartElement("format");
            	this.mWriter.writeAttribute("ID",line.getID());
            	this.mWriter.writeAttribute("type",line.getType().name());
            	if(line.isFixedCount()) this.mWriter.writeAttribute("count",String.valueOf(line.getCount()));
            	this.mWriter.writeCharacters(line.getDescription());
            	this.mWriter.writeEndElement();
        	  	}
        	this.mWriter.writeEndElement();
        	
        	//FILTER
        	this.mWriter.writeStartElement("filters-list");
        	for ( VCFFilterHeaderLine line : header.getFilterLines() )
        	  	{
            	this.mWriter.writeStartElement("filter");
            	this.mWriter.writeAttribute("ID",line.getID());
            	this.mWriter.writeCharacters(line.getValue());
            	this.mWriter.writeEndElement();
        	  	}
        	this.mWriter.writeEndElement();

        	//CONTIGS
        	this.mWriter.writeStartElement("contigs-list");
        	for ( VCFContigHeaderLine line : header.getContigLines() )
        	  	{
            	this.mWriter.writeStartElement("contig");
            	this.mWriter.writeAttribute("ID",line.getID());
            	this.mWriter.writeAttribute("index",String.valueOf(line.getContigIndex()));
            	this.mWriter.writeEndElement();
        	  	}
        	this.mWriter.writeEndElement();
        	
        	//SAMPLES
        	this.mWriter.writeStartElement("samples-list");
        	for (int i=0;i< header.getSampleNamesInOrder().size();++i )
        	  	{
            	this.mWriter.writeStartElement("sample");
            	this.mWriter.writeAttribute("id",String.valueOf(i+1));
            	this.mWriter.writeCharacters(header.getSampleNamesInOrder().get(i));
            	this.mWriter.writeEndElement();
        	  	}
        	this.mWriter.writeEndElement();

        	this.mWriter.writeEndElement();//head
        	this.mWriter.writeStartElement("body");
        	this.mWriter.writeStartElement("variations");
        }
        catch (XMLStreamException e)
    		{
        	throw new ReviewedStingException("IOException writing the VCF/XML header to " + super.getStreamName(), e);
    		}

    	}
	
	@Override
    public void add(VariantContext vc)
    	{	
        try
	        {
        	super.add(vc);

            Map<Allele, String> alleleMap = buildAlleleMap(vc);
             
	        this.mWriter.writeStartElement("variation");
	        
	        this.mWriter.writeAttribute("chrom",vc.getChr());
	        this.mWriter.writeAttribute("pos",String.valueOf(vc.getStart()));

	        
            String ID = vc.getID();
            if(!(ID==null || ID.isEmpty() || ID.equals(".")))
            	{
    	        this.mWriter.writeStartElement("id");
    	        this.mWriter.writeCharacters(ID);
    	        this.mWriter.writeEndElement();//id
            	}
            
	        this.mWriter.writeStartElement("id");
	        this.mWriter.writeCharacters(ID);
	        this.mWriter.writeEndElement();//body

	        this.mWriter.writeStartElement("ref");
	        this.mWriter.writeCharacters( vc.getReference().getDisplayString());
	        this.mWriter.writeEndElement();


	        if ( vc.isVariant() )
	        	{
                for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
                	{
                    Allele altAllele = vc.getAlternateAllele(i);
        	        this.mWriter.writeStartElement("alt");
        	        this.mWriter.writeCharacters(altAllele.getDisplayString());
        	        this.mWriter.writeEndElement();
                	}
	        	}
	        
	        
	        this.mWriter.writeEndElement();//variation
	        }
	    catch(XMLStreamException err)
	    	{
	    	throw new ReviewedStingException("Cannot close XMLStream",err);
	    	}
		}

	
	@Override
    public void close()
    	{
        super.close();
        try
	        {
	        this.mWriter.writeEndElement();//variations
	        this.mWriter.writeEndElement();//body
	        this.mWriter.writeEndElement();//vcf
	        this.mWriter.flush();
	        this.mWriter.close();
	        }
        catch(XMLStreamException err)
        	{
        	throw new ReviewedStingException("Cannot close XMLStream",err);
        	}
        try
	        {
	        getOutputStream().close();
	        }
	      catch(Throwable err)
        	{
        	throw new ReviewedStingException("Cannot close ouputstream",err);
        	}
    	}
	}