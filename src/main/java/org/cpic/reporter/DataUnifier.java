package org.cpic.reporter;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.cpic.reporter.model.CPICException;
import org.cpic.reporter.model.CPICinteraction;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.DiplotypeCall;
import org.cpic.reporter.resultsJSON.Gene;

public class DataUnifier {
    List<DiplotypeCall> calls;
    Map<String, List<CPICException>> exceptions;
    List<CPICinteraction> drugGenes;
    
   
    
   public DataUnifier( List<DiplotypeCall> calls, 
                                   Map<String, List<CPICException>> matches,
                                   List<CPICinteraction> drugGenes){
        this.calls = calls;
        this.exceptions = matches;
        this.drugGenes = drugGenes;
    }
    
    public void findMatches(){
        
        List<Gene> callSetToReturn = new ArrayList<Gene>();
        ExceptionMatcher matchTest = new ExceptionMatcher();
        
        for(DiplotypeCall call : calls){
            
            Gene gene = new Gene(call);
            
            if( exceptions.containsKey(call.getGene()) ){
                for( CPICException exception : exceptions.get(call.getGene() ) ){
                   if( matchTest( gene, exception.getRule_name()) ){
                       gene.addException(exception);
                   }
                       
                }
                
            }
        }
    }

}
