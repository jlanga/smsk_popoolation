import java.util.regex.Pattern;
import java.util.*;



public class MpileupParser {
	private final int minimumQuality;
	private final QualityEncoding qualityEncoding;
	
	public MpileupParser(QualityEncoding qualityEncoding, int minimumQuality)
	{
		this.minimumQuality=minimumQuality;
		this.qualityEncoding=qualityEncoding;
	}
	
	/*
	 * Parse a line from a synchronized line
	 */
	public  SyncLine parseLine(String line)
	{
		String[] c=line.split("\t");
		String chromosome=c[0];
		int position=Integer.parseInt(c[1]);
		char refChar=c[2].charAt(0);
		
		ArrayList<SyncPop> populations=new ArrayList<SyncPop>();
		for(int i=3; i<c.length; i+=3)
		{
			// c[i] == coverage
			int cov=Integer.parseInt(c[i]);
			SyncPop sp;
			if(cov > 0)
			{
				String seq=c[i+1];
				String qual=c[i+2];
				sp =getSyncPop(seq,qual,refChar);
			}
			else
			{
				sp=new SyncPop(0,0,0,0,0,0);
			}

			populations.add(sp);
		}
		return new SyncLine(chromosome,position,refChar,populations);
	}
	
	private SyncPop getSyncPop(String seq, String qual,char refchar)
	{
		if(seq.equals('-')) return new SyncPop(0,0,0,0,0,0);
		char[] seqca=seq.toCharArray();
		char[] qualca = qual.toCharArray();
		
		ArrayList<Character> purgedseq=purgeSequence(seqca);
		assert(purgedseq.size()==qualca.length);
		
		
		int countA=0; int countT=0; int countC=0; int countG=0; int countN=0; int countDel=0;
		for(int i=0; i<purgedseq.size(); i++)
		{
			char s=purgedseq.get(i);
			char q=qualca[i];
			int bqual=translateQuality(q);
			
			// Discard bases not having the minimum quality
			if(bqual < this.minimumQuality) continue;
			
			if(s== '.' || s==',') s=refchar;
			
			if(		s=='A' || s=='a') countA++;
			else if(s=='T' || s=='t') countT++;
			else if(s=='C' || s=='c') countC++;
			else if(s=='G' || s=='g') countG++;
			else if(s=='N' || s=='n') countN++;
			else if(s=='*') 			countDel++;
			else throw new IllegalArgumentException("Do not recognise character: " +s);
		}
		
		return new SyncPop(countA,countT,countC,countG,countN,countDel);
	}
	
	
	
	private int translateQuality(char c)
	{
		int raw=(int)c;
		int qual=0;
		if(this.qualityEncoding==QualityEncoding.Sanger)
		{
			qual= (raw-33);
		}
		else if(this.qualityEncoding==QualityEncoding.Illumina)
		{
			qual= (raw-64);
		}
		else
		{
			throw new IllegalArgumentException("Unknown quality encoding");
		}
		
		if(qual<0) throw new IllegalArgumentException("Found quality lower than zero; This is impossible; Please use the proper quality encoding (sanger|illumina)");
		return qual;
		
	}
	
	/*
	 * Center piece of parsing
	 * Get rid of all that 'schmarn' in the sequence of the mpileup
	 */
	private ArrayList<Character> purgeSequence(char[] toPurge)
	{
		ArrayList<Character> purged=new ArrayList<Character>();
		
		// # get rid of the crap present in the sequence line
	    // $nucs=~s/[-+](\d+)(??{"[ACGTNacgtn]{$1}"})//g; # I offer a beer for anyone who understand this line, I am a genius!!!!
	    // $nucs=~s/\^.//g;
	    //$nucs=~s/\$//g;
	    //$nucs=~s/[.]/uc($rc)/eg;
	    //$nucs=~s/[,]/lc($rc)/eg;
		int i =0;
		while(i<toPurge.length)
		{
			// active character= ac
			char ac=toPurge[i];
			if(ac == '+' || ac=='-')
			{
				ArrayList<Character> digits=new ArrayList<Character>();
				
				// Find digits
				int k=1;
				while(Character.isDigit(toPurge[i+k]))
				{
					digits.add(toPurge[i+k]);
					k++;
				}
				
				StringBuilder sb=new StringBuilder();
				for(char c: digits) sb.append(c);
				int indelsize=Integer.parseInt(sb.toString());
				
				//s/[-+](\d+)(??{"[ACGTNacgtn]{$1}"})//g; # I offer a beer for anyone who understand this line, I am a genius!!!!
				// [+-]  = 1
				// (\d+) = digits.size()
				// [ACGTNacgtn]{$1} = indelsize
				int ignoresize=1 + digits.size()+ indelsize;
				i+=ignoresize;
				
			}
			else if(ac=='^')
			{
				i+=2;
			}
			else if(ac=='$')
			{
				i++;
			}
			else
			{
				purged.add(ac);
				i++;
			}
		}
		return purged;
	}
}
