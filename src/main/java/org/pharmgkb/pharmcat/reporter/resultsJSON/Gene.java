package org.pharmgkb.pharmcat.reporter.resultsJSON;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.json.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.CPICException;


public class Gene {

    private String gene;
    Set<String> diplotypes;
    List<CPICException> exceptList = new ArrayList();
    List<Interaction> drugGene = new ArrayList();;

    public Gene(GeneCall call) {
        this.gene = call.getGene();
      diplotypes = call.getDiplotypes().stream()
          .map(DiplotypeMatch::getName)
          .collect(Collectors.toSet());
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

