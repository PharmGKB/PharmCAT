package org.pharmgkb.pharmcat.reporter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.pharmgkb.pharmcat.haplotype.model.json.GeneCall;
import org.pharmgkb.pharmcat.reporter.io.JsonFileLoader;
import org.pharmgkb.pharmcat.reporter.model.CPICException;
import org.pharmgkb.pharmcat.reporter.model.CPICinteraction;
import org.pharmgkb.pharmcat.reporter.resultsJSON.Gene;
import org.pharmgkb.pharmcat.reporter.resultsJSON.MultiGeneInteraction;

//import org.slf4j.Logger;
//import org.slf4j.LoggerFactory;

/**
 * This is contains the main function for running the CPIC reporting tool.
 * 
 * As of today  (4-23-2016) this tool is still in development currently this tool 
 * there are many broken, missing, and commented out pieces of this code due to
 * ongoing work to produce a working output.  
 * 
 * Please contact me (Greyson Twist) on slack or at gtwist@cmh.edu if there are 
 * questions and I will try to improve documentation to make thing more clear
 * 
 * 
 * @author greytwist
 *
 */


public class Reporter {


    /**
     * Logging instance
     */
   //private Logger logger = LoggerFactory.getLogger( Reporter.class );


	String locationOfResources = "C:\\Users\\greytwist\\Desktop\\CLONES\\PharmCAT\\";
   /**
    * Exception list formated as json
    */
   String exceptionPath = "resources\\cpic_exceptions\\exceptions.json"; //TODO don't do this, done for wiring purposes only
   private File exception = new File(locationOfResources + exceptionPath);

   /**
    * Drug Gene interaction json
    */

   //TODO don't do any of this load from props file like a read engineer
   String test1 = "resources\\dosing_guidelines\\CPIC_Guideline_for_citalopram_escitalopram_and_CYP2C19.json";
   String test2 = "resources\\dosing_guidelines\\CPIC_Guideline_for_clopidogrel_and_CYP2C19.json";
   String test3 = "resources\\dosing_guidelines\\CPIC_Guideline_for_sertraline_and_CYP2C19.json";
   private List<File> interactions = new ArrayList<File>();

    /**
     * Configuration properties loaded from file
     */
   // private Properties props;

    /**
     * File in
     * TODO CLEAN THIS UP FOR TEST BUILDING AND WIRING ONLY
     */
    String multiGeneFile = "resources\\json_out_example\\multi_gene.json";
    private File inFile = new File(locationOfResources + multiGeneFile);

    /**
     * File to write final results to
     */
   // private File outFile;



    //parse command line options
    private Reporter( CommandLine cmdline )  throws IOException {
        interactions.add( new File(locationOfResources + test1));
        interactions.add( new File(locationOfResources + test2));
        interactions.add( new File(locationOfResources + test3));

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

        /*if( args.length == 0 ) {
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

    	List<Gene> geneListToReport = new ArrayList<Gene>();
    	List<MultiGeneInteraction> multiInteractionsToReport = new ArrayList<MultiGeneInteraction>(); 
        //loadRequiredFiles(); TODO undelete this and use actual args and real code for plumbing the system

    	//Generate class used for loading JSON into 
        JsonFileLoader loader = new JsonFileLoader();

        //Load the haplotype json, this is pointed at a test json and will likely break when meeting real 
        // requiring some if not all rewriting
        List<GeneCall> calls = loader.loadHaplotypeGeneCalls(this.inFile);
        
        //Load the gene rule exception json
        Map<String, List<CPICException>> exceptions = loader.loadExceptions(this.exception);
        
        //Load the gene drug interaction list. This currently only handles single gene drug interactions and will require updating to handle multi gene drug interaction
        List<CPICinteraction> drugGenes = loader.loadDrugGeneRecommendations(this.interactions);

        //This is the primary work flow for generating the report where calls are matched to exceptions and drug gene interactions based on reported haplotypes
        DataUnifier checker = new DataUnifier(calls, exceptions, drugGenes); // prime with data
        checker.findMatches(geneListToReport, multiInteractionsToReport ); // run the actual comparison



         //TODO print results here!!!!!

       // logger.info( "Complete" );
    }

    /*private void loadRequiredFiles(){
        String exceptionLoc = props.getProperty("CPIC.reporter.exception");
        this.exception = new File( exceptionLoc );

        String interactionLoc = props.getProperty("CPIC.reporter.guidlines");
        this.interactions = new File(interactionLoc);
    }*/






}
