package org.cpic.reporter.resultsJSON;

import java.util.List;
import java.util.Set;

import org.cpic.reporter.model.CPICException;
import org.cpic.reporter.model.CPICinteraction;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.DiplotypeCall;

public class Gene {
    
    private String gene;
    Set<String> diplotypes;
    List<CPICException> exceptlist; 
    List<Interaction> drugGene;
    
    public Gene( DiplotypeCall call){
        this.gene = call.getGene();
        this.diplotypes = call.getDiplotypes();
    }
    
    

}
