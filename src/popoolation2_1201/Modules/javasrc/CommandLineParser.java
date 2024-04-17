import java.util.Arrays;
import java.util.LinkedList;


public class CommandLineParser 
{
	
	public static void parseArguments(String[] arguments,java.util.logging.Logger logger)
	{
	

    LinkedList<String> args=new LinkedList<String>(Arrays.asList(arguments));
    String input="";
    String output="";
    QualityEncoding qualEncoding=QualityEncoding.Sanger;
    int minimumQual=20;
    int threads=1;
    int blockSize=50000;
    boolean help=false;

    
    // Parse the command line arguments
    // order does not matter
    while(args.size()>0)
    {
        String cu=args.remove(0);
        if(cu.equals("--input"))
        {
            input=args.remove(0);
        }
        else if(cu.equals("--output"))
        {
            output=args.remove(0);
        }        
        else if(cu.equals("--fastq-type"))
        {
            String rawQualEncoding=args.remove(0);
            qualEncoding=translateQualityEncoding(rawQualEncoding);
        }
        else if(cu.equals("--min-qual"))
        {
            minimumQual=Integer.parseInt(args.remove(0));
        }
        else if(cu.equals("--threads"))
        {
        	threads=Integer.parseInt(args.remove(0));
        }
        else if(cu.equals("--help"))
        {
        	help=true;
        }
        else
        {
        	throw new IllegalArgumentException("Do not recognise argument "+cu);
        }
    }
        
        // If help was requested -> show help message and exit
       if(help)
       {
        	System.out.print(getHelpMessage());
        	System.exit(1);
       }
        
        
        Mpileup2SyncFramework frame= new Mpileup2SyncFramework(input,output,qualEncoding,minimumQual,blockSize,threads,logger);
        frame.run();
        
        
    }

	
	private static QualityEncoding translateQualityEncoding(String rawQuality)
	{
		if(rawQuality.toLowerCase().equals("sanger")){
			return QualityEncoding.Sanger;
		}
		else if(rawQuality.toLowerCase().equals("illumina"))
		{
			return QualityEncoding.Illumina;
		}
		else
		{
			throw new IllegalArgumentException("Do not recognise quality encoding "+ rawQuality);
		}
	}
	
	private static String getHelpMessage(){
		StringBuilder sb=new StringBuilder();
		sb.append("mpileup2sync: converts a mpileup file into a synchronized file\n");
		sb.append("--input					the input file in the mpileup format; mandatory\n");
		sb.append("--output					the output file, will be a synchronized file; mandatory\n");
		sb.append("--fastq-type			the encoding of the quality sequences; sanger|illumina; default=sanger\n");
		sb.append("--min-qual			the minimum base quality to filter for; default=0\n");
		sb.append("--threads					number of threads to use; default=1\n");
		sb.append("--help					display the help pages\n");
		return sb.toString();
	}
	
}
