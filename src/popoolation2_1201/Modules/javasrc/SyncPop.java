
public class SyncPop {
	//A-count:T-count:C-count:G-count:N-count:deletion-count
	private final int countA;
	private final int countT;
	private final int countC;
	private final int countG;
	private final int countN;
	private final int countDel;
	public SyncPop(int countA, int countT, int countC, int countG, int countN, int countDel)
	{
		this.countA=countA;
		this.countT=countT;
		this.countC=countC;
		this.countG=countG;
		this.countN=countN;
		this.countDel=countDel;
	}
	
	
	public int countA()
	{
		return this.countA;
	}
	public int countT()
	{
		return this.countT;
	}
	public int countC()
	{
		return this.countC;
	}
	public int countG()
	{
		return this.countG;
	}
	public int countN()
	{
		return this.countN;
	}
	public int countDel()
	{
		return this.countDel;
	}

}
