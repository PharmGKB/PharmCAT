package org.pharmgkb.pharmcat.reporter.resultsJSON;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
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
