package org.cpic.reporter.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import com.google.gson.Gson;
import org.cpic.reporter.model.CPICException;
import org.cpic.reporter.model.CPICExceptionList;
import org.cpic.reporter.model.CPICinteraction;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.DiplotypeCall;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.HaplotyperResult;
public class JsonFileLoader {

    Gson gson = new Gson();

    public List<DiplotypeCall> loadHaplotypeGeneCalls( File haplotypeCalledFile ) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader( haplotypeCalledFile ));
        HaplotyperResult calls = gson.fromJson(br, HaplotyperResult.class);
        br.close();
        return calls.getDiplotypeCalls();

    }
    
    /**
     * This turns the json exception list into a fragmented list under the header
     * of the relevant gene.
     * FIXME dumb assumption that the exception are for single gene issues and not for 
     * multi gene interactions
     */
    public Map<String, List<CPICException>> loadExceptions( File exceptions )throws IOException {
        Map<String, List<CPICException>> matcher = new HashMap<String, List<CPICException>>();
        BufferedReader br = new BufferedReader(new FileReader( exceptions ));
        CPICExceptionList except = gson.fromJson(br, CPICExceptionList.class);
        for (CPICException rule : except.getRules()) {
            if( matcher.containsKey(rule.getGene())){
                matcher.get(rule.getGene()).add(rule);
            } else {
                List<CPICException> parts = new ArrayList<CPICException>();
                parts.add(rule);
                matcher.put(rule.getGene(), parts);
            }
        }

        br.close();

        return matcher;

    }

    /**
     * DANGER WILL ROBINSON! A DUMB ASSUMPTION HAS BEEN MADE!
     * FIXME : assumption of single gene and single drug interaction here
     *
     * Multi gene drugs need alteration in this file to work but this will be 
     * acceptable for hackathon standards but should be fixed as part of a the 
     * next version.  This is a critical flaw in the design
     *
     */
    public Map<String, List<CPICinteraction>> loadDrugGeneRecommendations( List<File> interactions ) throws IOException {
        Map<String,List<CPICinteraction>> drugGenes = new HashMap<String, List<CPICinteraction>>();
        for( File interact: interactions ){
           
            BufferedReader br = new BufferedReader(new FileReader( interact ));
            CPICinteraction act = gson.fromJson(br, CPICinteraction.class);
            String geneSymbol = act.getRelatedGenes().get(0).getSymbol(); //FIXME setting index to 0 will only allow the first gene if there are multiple
            
            if( drugGenes.containsKey(geneSymbol)){
                drugGenes.get(geneSymbol).add(act);
            } else {
                List<CPICinteraction> parts = new ArrayList<CPICinteraction>();
                parts.add(act);
                drugGenes.put(geneSymbol, parts);
            }
            br.close();
        }
        
        return drugGenes;

    }
}
