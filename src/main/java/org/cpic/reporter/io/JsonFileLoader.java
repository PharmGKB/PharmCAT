package org.cpic.reporter.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import org.cpic.reporter.model.CPICinteraction;
import org.cpic.reporter.model.CpicException;
import org.cpic.reporter.model.CpicExceptionList;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.DiplotypeCall;
import org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON.HaplotypeCallsMultiGene;

import com.google.gson.Gson;

public class JsonFileLoader {
    
    Gson gson = new Gson();

    public void loadHaplotypeGeneCalls( File haplotypeCalledFile ) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader( haplotypeCalledFile ));
        HaplotypeCallsMultiGene calls = gson.fromJson(br, HaplotypeCallsMultiGene.class);
        System.out.println( "Haplotype Caller test" );
        for (DiplotypeCall call : calls.getDiplotypeCalls() ) {
            System.out.println(call.getGene());
        }
        br.close();
        
    }
    
    public void loadExceptions( File exceptions )throws IOException {
        BufferedReader br = new BufferedReader(new FileReader( exceptions ));
        CpicExceptionList except = gson.fromJson(br, CpicExceptionList.class);
        System.out.println( "Exception test" );
        for (CpicException rule : except.getRules()) {
            System.out.println(rule.getGene());
          
        }

        br.close();
        
    }
    
    public void loadDrugGeneRecommendations( List<File> interactions ) throws IOException {
        for( File interact: interactions ){
            BufferedReader br = new BufferedReader(new FileReader( interact ));
            CPICinteraction act = gson.fromJson(br, CPICinteraction.class);
            System.out.println("Interactions List");
            System.out.println( act.getRelatedGenes().get(0).getSymbol() );
            br.close();
        }
        
    }
}
