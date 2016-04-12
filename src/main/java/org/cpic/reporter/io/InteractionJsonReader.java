package org.cpic.reporter.io;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;

import org.cpic.reporter.model.CPICinteraction;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;


public class InteractionJsonReader {
    
    public void load( File fileToRead ){
        try(Reader reader = new InputStreamReader(CPICinteraction.class.getResourceAsStream( fileToRead ), "UTF-8")){
            Gson gson = new GsonBuilder().create();
            CPICinteraction p = gson.fromJson(reader, CPICinteraction.class);
            System.out.println(p);
        }

}
