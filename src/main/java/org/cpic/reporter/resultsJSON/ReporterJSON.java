package org.cpic.reporter.resultsJSON;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import com.google.common.base.Preconditions;
import com.google.gson.FieldNamingPolicy;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

public class ReporterJSON {
	
	public static void main(String[] args) {
		JsonObject jsonObject = new JsonObject();
//		Gson gson = new Gson();\
		
		ReporterResult results = new ReporterResult();
		Gson gson = new GsonBuilder().setPrettyPrinting().create();
				
		
		//JsonElement jsonElement = ;
		jsonObject.add("gene", gson.toJsonTree(ReporterResult.class.getName())); 
//		jsonObject.add("diplotypes", gson.toJsonTree(ReporterResult.class); 
    }
}
