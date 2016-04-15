package org.cpic.reporter.resultsJSON;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import org.cpic.haplotype.model.json.DiplotypeCall;
import org.cpic.reporter.model.CPICException;

public class Gene {

    private String gene;
    Set<String> diplotypes ;
    List<CPICException> exceptList = new ArrayList();
    List<Interaction> drugGene = new ArrayList();;

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
}

