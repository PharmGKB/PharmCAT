package org.pharmgkb.pharmcat.reporter.resultsJSON;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import org.pharmgkb.pharmcat.haplotype.model.json.DiplotypeCall;
import org.pharmgkb.pharmcat.reporter.model.CPICException;


public class Gene {

    private String gene;
    Set<String> diplotypes ;
    List<CPICException> exceptList = new ArrayList<CPICException>();
    List<Interaction> drugGene = new ArrayList<Interaction>();;

    public Gene( DiplotypeCall call){
        this.gene = call.getGene();
        this.diplotypes = call.getDiplotypes();
        
        
    }

    public void addException( CPICException except ){
        exceptList.add(except);
    }

    public Set<String> getDips(){
        return diplotypes;
    }

    public String getGene(){
        return gene;
    }

    public List<CPICException> getExceptionList(){
        return exceptList;
    }

    public void addInteraction( Interaction interact) {
        drugGene.add(interact);
    }
    
    public List<Interaction> getInteractionList(){
    	return drugGene;
    }
}

