import java.util.logging.Level;


public class Mpileup2syncMain {
	
    public static void main(String[] args) {
    	
        java.util.logging.Logger logger=java.util.logging.Logger.getLogger("Gowinda Logger");
        java.util.logging.ConsoleHandler gowhandler =new java.util.logging.ConsoleHandler();
        gowhandler.setLevel(Level.INFO);
        gowhandler.setFormatter(new LogFormatter());
        logger.addHandler(gowhandler);
        logger.setUseParentHandlers(false);
        logger.setLevel(Level.ALL);

    	
    	CommandLineParser.parseArguments(args,logger);
    	System.exit(1);
    
    }

}
