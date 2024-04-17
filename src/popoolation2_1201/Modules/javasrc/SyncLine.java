
import java.util.*;

public class SyncLine {
	private final String chromosome;
	private final int position;
	private final char refChar;
	private final ArrayList<SyncPop> populations;

	public SyncLine(String chromosome, int position, char refChar, ArrayList<SyncPop> populations)
	{
		this.chromosome=chromosome;
		this.position=position;
		this.refChar=refChar;
		this.populations=new ArrayList<SyncPop>(populations);
	}
	
	public String chromosome()
	{
		return this.chromosome;
	}
	public int position()
	{
		return this.position;
	}
	public char referenceCharacter()
	{
		return this.refChar;
	}
	public ArrayList<SyncPop> populations()
	{
		return new ArrayList<SyncPop>(this.populations);
	}
}

