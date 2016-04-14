package org.cpic.reporter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.cpic.reporter.io.JsonFileLoader;
//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;


public class Reporter {


    /**
     * Logging instance
     */
   //private Logger logger = LoggerFactory.getLogger( Reporter.class );


   /**
    * Exception list formated as json
    */
   String exceptionPath = "/gpfs/data/home/gtwist/tmp/CPIC/cpic-annotator/resources/cpic_exceptions/exceptions.json"; //TODO don't do this, done for wiring purposes only
   private File exception = new File(exceptionPath);
   
   /**
    * Drug Gene interaction json
    */
   
   //TODO don't do any of this load from props file like a read engineer
   String test1 = "/gpfs/data/home/gtwist/tmp/CPIC/cpic-annotator/resources/CPIC_Guideline_for_citalopram_escitalopram_and_CYP2C19.json";
   String test2 = "/gpfs/data/home/gtwist/tmp/CPIC/cpic-annotator/resources/CPIC_Guideline_for_clopidogrel_and_CYP2C19.json";
   String test3 = "/gpfs/data/home/gtwist/tmp/CPIC/cpic-annotator/resources/CPIC_Guideline_for_sertraline_and_CYP2C19.json";
   private List<File> interactions = new ArrayList<File>();
   
    /**
     * Configuration properties loaded from file
     */
   // private Properties props;

    /**
     * File in
     * TODO CLEAN THIS UP FOR TEST BUILDING AND WIRING ONLY
     */
    String multiGeneFile = "/gpfs/data/home/gtwist/tmp/CPIC/cpic-annotator/resources/cpic_exceptions/exceptions.json";
    private File inFile = new File(multiGeneFile);

    /**
     * File to write final results to
     */
   // private File outFile;



    //parse command line options
    private Reporter( CommandLine cmdline )  throws IOException {
        interactions.add( new File(test1));
        interactions.add( new File(test2));
        interactions.add( new File(test3));
        
        //File propsFile = new File( cmdline.getOptionValue( "conf" ) );
       // props = readConfFile( propsFile );

        //String out_location = cmdline.getOptionValue( "outFile" );
        //this.outFile = new File( out_location );

       // String in_file = cmdline.getOptionValue("inFile");
        //this.inFile = new File( in_file );


    }

    /**
     * Read config file into props object
     */
    /*private Properties readConfFile( File file ) throws IOException {
        Properties props = new Properties();
        props.load( FileUtils.openInputStream( file ) );
        return props;
    }*/

    public static Options createCommandLineOptions() {
        Options options = new Options();

        Option infile = new Option("inFile", true, "required - haplotype result to generate report");
        options.addOption(infile);

        Option conf = new Option( "conf", true, "required - configuration file" );
        options.addOption( conf );


        Option outfile = new Option( "outFile", true, "required - location to write final file" );
        options.addOption( outfile );

        return options;
    }

    public static void main(String[] args) throws Exception {
        Options options = createCommandLineOptions();
        CommandLine cmdline = new DefaultParser().parse( options,  args );

        /*/8if( args.length == 0 ) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( "CPIC-Reporter", options );
            System.exit( 0 );
        } else if ( !cmdline.hasOption("conf") ||
                    !cmdline.hasOption( "in" ) ||
                    !cmdline.hasOption("out") ) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp( "CPIC-Reporter", options );
            System.out.println( "Error in options in, out, conf" +
                    " are required variables" );
            System.exit( 1 );

        } else if ( !cmdline.hasOption("inputBam")){
            System.out.println( "WARNING: No bam file provided" );
        }*/

        
        //if minimal required parameters are set parse command line
        Reporter report = new Reporter( cmdline );
        //run genotyper workflow
        report.run();
    }

    private void run() throws IOException {
        
        //loadRequiredFiles(); TODO undelete this and use actual args and real code for plumbing the system
        
        JsonFileLoader loader = new JsonFileLoader();
        loader.loadHaplotypeGeneCalls(this.inFile);
        loader.loadExceptions(this.exception);
        loader.loadDrugGeneRecommendations(this.interactions);
        
        


       // logger.info( "Complete" );
    }
    
    /*private void loadRequiredFiles(){
        String exceptionLoc = props.getProperty("CPIC.reporter.exception");
        this.exception = new File( exceptionLoc );
        
        String interactionLoc = props.getProperty("CPIC.reporter.guidlines");
        this.interactions = new File(interactionLoc);
    }*/
    





}
