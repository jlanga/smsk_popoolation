import java.io.*;
public class BatchWriter {
	private  BufferedWriter bw;
	private final String outputFile;
	
	
	public BatchWriter(String outputFile)
	{
	
		try
		{
			this.bw = new BufferedWriter(new FileWriter(outputFile));
		}
		catch(IOException e)
		{
			e.printStackTrace();
			System.exit(0);
		}

		this.outputFile=outputFile;
	}
	
	public void writeBatch(String[] batch)
	{
		try{
			for(String s: batch)
			{
				bw.write(s+"\n");
			}			
		}
		catch(IOException e)
		{
			e.printStackTrace();
			System.exit(0);
		}

	}
	
	public void close()
	{
		try
		{
			bw.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
			System.exit(0);
		}
	}

}
