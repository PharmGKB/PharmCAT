package org.cpic.reporter;

import java.util.List;
import java.util.Map;

import org.cpic.reporter.model.CPICException;
import org.cpic.reporter.model.CPICinteraction;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.DiplotypeCall;

public class RecommendationMatcher {
    List<DiplotypeCall> calls;
    Map<String, CPICException> matches;
    Map<String, CPICinteraction> drugGenes;
    
    public void RecoommendationMatcher( List<DiplotypeCall> calls, 
                                   Map<String, CPICException> matches,
                                   Map<String, CPICinteraction> drugGenes){
        this.calls = calls;
        this.matches = matches;
        this.drugGenes = drugGenes;
    }

}
