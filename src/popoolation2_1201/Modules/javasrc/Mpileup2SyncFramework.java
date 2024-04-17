import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.util.logging.Logger;

public class Mpileup2SyncFramework {
	private final String input;
	private final String output;
	private final QualityEncoding qualityEncoding;
	private final int minimumQuality;
	private final int blockSize;
	private final int threads;
	private final Logger logger;
	
	
	public Mpileup2SyncFramework(String input, String output, QualityEncoding qualityEncoding, int minimumQuality, int blockSize, int threads, Logger logger)
	{
        if(!new File(input).exists()){throw new IllegalArgumentException("Input file does not exist "+ input);}
        if(threads<1) throw new IllegalArgumentException("Number of threads needs to be larger than zero");
        if(blockSize<1) throw new IllegalArgumentException("Block size must be larger than 1; Should be >>1");
        if(minimumQuality<0) throw new IllegalArgumentException("Minimum quality must be larger than zero");
        
        try
        {
        	new FileWriter(output);
        }
        catch(IOException e)
        {
        	throw new IllegalArgumentException("Can not create output file:" +output);
        }
		
		this.input=input;
		this.output=output;
		this.qualityEncoding=qualityEncoding;
		this.minimumQuality=minimumQuality;
		this.blockSize=blockSize;
		this.threads=threads;
		this.logger=logger;
	}
	
	
	public void run()
	{
		BatchReader br=new BatchReader(this.input,this.blockSize);
		BatchWriter bw=new BatchWriter(this.output);
		ExecutorService executor=Executors.newFixedThreadPool(threads);
		MpileupParser parser =new MpileupParser(this.qualityEncoding,this.minimumQuality);
		
		String[] batch;
		while((batch=br.getBatch()).length>0)
		{
			String[] syncBatch=new BatchProcessor(batch,parser,executor).processBatch();
			assert(syncBatch.length==batch.length);
			bw.writeBatch(syncBatch);
		}
		
		br.close();
		bw.close();
		
	}

}
