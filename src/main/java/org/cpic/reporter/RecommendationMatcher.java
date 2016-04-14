package org.cpic.reporter;

import java.util.List;
import java.util.Map;

import org.cpic.reporter.model.CPICException;
import org.cpic.reporter.model.CPICinteraction;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.DiplotypeCall;

public class RecommendationMatcher {
    List<DiplotypeCall> calls;
    Map<String, List<CPICException>> exceptions;
    Map<String, CPICinteraction> drugGenes;
    
    public void RecoommendationMatcher( List<DiplotypeCall> calls, 
                                   Map<String, List<CPICException>> matches,
                                   Map<String, CPICinteraction> drugGenes){
        this.calls = calls;
        this.exceptions = matches;
        this.drugGenes = drugGenes;
    }
    
    public void findMatches(){
        for(DiplotypeCall call : calls){
            System.out.println("Gene: " + call.getGene());
            
            if( exceptions.containsKey(call.getGene()) ){
                
            } else {
                
            }
        }
    }

}
