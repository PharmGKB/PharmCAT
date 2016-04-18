package org.pharmgkb.pharmcat.reporter.development;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import com.google.gson.Gson;
import org.pharmgkb.pharmcat.reporter.model.CalledType;


public class TestJsonParserForHap {

    public static void main(String[] args) throws Exception {
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader("/gpfs/data/home/gtwist/tmp/CPIC/cpic-annotator/resources/json_out_example/CYP2C19_multiple.json"));
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        Gson gson = new Gson();
        CalledType call = gson.fromJson(br, CalledType.class);
        System.out.println( call.getMetaData().getCpicAnnotatorBuild() );
      
                
    }
          
    
    

}
