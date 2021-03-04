import java.io.*;
import java.util.*;

public class BatchReader {
	private String inputFile;
	private int blockSize;
	private BufferedReader sr;
	public BatchReader(String inputFile, int blockSize)
	{
		try{
			sr=new BufferedReader(new FileReader(inputFile));
		}
		catch(FileNotFoundException e)
		{
			e.printStackTrace();
			System.exit(0);
		}
		this.inputFile=inputFile;
		this.blockSize=blockSize;
	}
	
	
	public String[] getBatch()
	{
		ArrayList<String> toret=new ArrayList<String>();
		try
		{
			String line;
			int eCount=0;
			while(eCount < this.blockSize && (line=sr.readLine())!=null)
			{
				toret.add(line);
				eCount++;
			}		
		}
		catch(IOException e)
		{
			e.printStackTrace();
			System.exit(0);
		}
		
		
		// Convert to String[]
		String[] tr=new String[toret.size()];
		int counter=0;
		for(String l: toret)
		{
			tr[counter]=l;
			counter++;
		}
		
		return  tr;
	}
	
	public void close()
	{
		try
		{
			sr.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
			System.exit(0);
		}
	}

}
