package org.cpic.reporter;

import java.util.List;
import java.util.Map;

import org.cpic.reporter.model.CPICException;
import org.cpic.reporter.model.CPICinteraction;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.DiplotypeCall;

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
        ExceptionMatcher matchTest = new ExceptionMatcher();
        for(DiplotypeCall call : calls){
            if( exceptions.containsKey(call.getGene()) ){
                for( CPICException exception : exceptions.get(call.getGene() ) ){
                    
                }
                
            }
        }
    }

}
