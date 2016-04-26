package org.pharmgkb.pharmcat.reporter.io;

import java.util.List;

import org.pharmgkb.pharmcat.reporter.resultsJSON.Gene;
import org.pharmgkb.pharmcat.reporter.resultsJSON.Interaction;
import org.pharmgkb.pharmcat.reporter.resultsJSON.MultiGeneInteraction;

public class ReporterWriter {

	public static void printResults(/*File outFile,*/ List<Gene> geneListToReport, List<MultiGeneInteraction> multiInteractionsToReport) {
		for( Gene geneToWrite : geneListToReport ){
			System.out.println( "Gene: " + geneToWrite.getGene() );
			System.out.println( "Called Diplotypes: " + geneToWrite.getGene() );
			
			int exceptCount = geneToWrite.getExceptionList().size();
			for( int i = 0; i <= exceptCount; i++ ){
				if(i == 0){
					System.out.println("Exceptions:");
				}
				System.out.println( "Rule: " + geneToWrite.getExceptionList().get(i).getRule_name() );
				System.out.println( "Version: " + geneToWrite.getExceptionList().get(i).getVersion() );
				System.out.println( "Matches: " + geneToWrite.getExceptionList().get(i).getMatches() );
				System.out.println( "Exception type: " + geneToWrite.getExceptionList().get(i).getException_type() );
				System.out.println( "Message: " + geneToWrite.getExceptionList().get(i).getMessage() );
				
			}
			
			int interactionCount = geneToWrite.getInteractionList().size();
			for( int i = 0; i <= interactionCount; i++ ){
				Interaction toWrite = geneToWrite.getInteractionList().get(i);
				System.out.println( "Name: " + toWrite.getName());
				System.out.println( "Source: " + toWrite.getSource());
				System.out.println( "SummaryHtml: " + toWrite.getSummaryHtml());
				System.out.println( "Html Text: " + toWrite.getTextHtml());
				
				
			}
			
		}
		
	}

}
