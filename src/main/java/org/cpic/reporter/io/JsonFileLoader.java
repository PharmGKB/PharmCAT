package org.cpic.reporter.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import org.cpic.reporter.model.CPICinteraction;
import org.cpic.reporter.model.CalledType;

import com.google.gson.Gson;

public class JsonFileLoader {

    public void loadHaplotypeGeneCalls( File haplotypeCalledFile ) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader( haplotypeCalledFile ));
        Gson gson = new Gson();
        CalledType call = gson.fromJson(br, CalledType.class);
        System.out.println( call.getMetaData().getCpicAnnotatorBuild() );
        br.close();
        
    }
    
    public void loadExceptions( File exceptions )throws IOException {
        BufferedReader br = new BufferedReader(new FileReader( exceptions ));
        Gson gson = new Gson();
       // CalledType call = gson.fromJson(br, CalledType.class);
        //System.out.println( call.getMetaData().getCpicAnnotatorBuild() );
        br.close();
        
    }
    
    public void loadDrugGeneRecommendations( List<File> interactions ) throws IOException {
        for( File interact: interactions ){
            BufferedReader br = new BufferedReader(new FileReader( interact ));
            Gson gson = new Gson();
            CPICinteraction act = gson.fromJson(br, CPICinteraction.class);
            System.out.println( act.getRelatedGenes().get(0).getName() );
            br.close();
        }
        
    }
}
