package org.cpic.reporter;

import java.io.File;
import java.io.IOException;
import java.util.Properties;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.io.FileUtils;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

public class Reporter {
    
    /**
     * Logging instance
     */
   private Log logger = LogFactory.getLog( Reporter.class );
    
    /**
     * Configuration properties loaded from file
     */
    private Properties props;
    
    /**
     * File in
     */
    private File inFile;
    
    /**
     * File to write final results to
     */
    private File outFile;
    

    
    //parse command line options
    private Reporter( CommandLine cmdline )  throws IOException {
        
        File propsFile = new File( cmdline.getOptionValue( "conf" ) );
        props = readConfFile( propsFile );
        
        String out_location = cmdline.getOptionValue( "outFile" );
        this.outFile = new File( out_location );
        
        String in_file = cmdline.getOptionValue("inFile");
        this.inFile = new File( in_file );
       
                                
    }
    
    /**
     * Read config file into props object
     */
    private Properties readConfFile( File file ) throws IOException {
        Properties props = new Properties();
        props.load( FileUtils.openInputStream( file ) );
        return props;
    }

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
        
        if( args.length == 0 ) {
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
        }
        
        //if minimal required parameters are set parse command line
        Reporter report = new Reporter( cmdline );
        //run genotyper workflow
        report.run();
    }
    
    private void run() throws IOException {
  
        
        logger.info( "Complete" );
    }

    
    
    
  

}
