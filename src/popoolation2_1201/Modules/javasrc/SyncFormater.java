
public class SyncFormater {
	
	public static String formatSync(SyncLine sl)
	{
		StringBuilder sb=new StringBuilder();
		sb.append(sl.chromosome()); 	sb.append("\t");
		sb.append(sl.position()); 		sb.append("\t");
		sb.append(sl.referenceCharacter());
		for(SyncPop sp : sl.populations())
		{
			sb.append("\t");
			String info= sp.countA()+":"+sp.countT()+":"+sp.countC()+":"+sp.countG()+":"+sp.countN()+":"+sp.countDel();
			sb.append(info);
		}
		return sb.toString();
	}

}
