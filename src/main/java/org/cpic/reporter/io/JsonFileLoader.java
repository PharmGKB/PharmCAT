package org.cpic.reporter.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.cpic.reporter.model.CPICException;
import org.cpic.reporter.model.CPICExceptionList;
import org.cpic.reporter.model.CPICinteraction;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.DiplotypeCall;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.HaplotypeCallsMultiGene;

import com.google.gson.Gson;
public class JsonFileLoader {
    
    Gson gson = new Gson();

    public List<DiplotypeCall> loadHaplotypeGeneCalls( File haplotypeCalledFile ) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader( haplotypeCalledFile ));
        HaplotypeCallsMultiGene calls = gson.fromJson(br, HaplotypeCallsMultiGene.class);
        br.close();
        return calls.getDiplotypeCalls();
        
    }
    
    public Map<String, CPICException> loadExceptions( File exceptions )throws IOException {
        Map<String, CPICException> matcher = new HashMap<String, CPICException>();
        BufferedReader br = new BufferedReader(new FileReader( exceptions ));
        CPICExceptionList except = gson.fromJson(br, CPICExceptionList.class);
        System.out.println( "Exception test" );
        for (CPICException rule : except.getRules()) {
           matcher.put(rule.getGene(), rule);
        }

        br.close();
        
        return matcher;
        
    }
    
    /**
     * DANGER WILL ROBINSON! A DUMB ASSUMPTION HAS BEEN MADE!
     * FIXME : assumption of single gene and single drug interaction here 
     * 
     * I (Greyson) did not see any actual data in any of the multi gene drug
     * guidelines however I made a note to ask about it tomorrow.  If this is 
     * still broken or you are the person that is responsible to fix this I apologize.
     * @return 
     * 
     */
    public Map<String, CPICinteraction> loadDrugGeneRecommendations( List<File> interactions ) throws IOException {
        Map<String, CPICinteraction> drugGenes = new HashMap<String, CPICinteraction>();
        for( File interact: interactions ){
            BufferedReader br = new BufferedReader(new FileReader( interact ));
            CPICinteraction act = gson.fromJson(br, CPICinteraction.class);
            drugGenes.put(act.getRelatedGenes().get(0).getSymbol(), act);
            br.close();
        }
        return drugGenes;
        
    }
}
