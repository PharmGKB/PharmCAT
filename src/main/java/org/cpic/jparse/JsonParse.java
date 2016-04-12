package src.main.java.org.cpic.jparse;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class JsonParse {
	public static void main(String[] args)  {
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader("/Users/averma/Downloads/dosingGuidelines.json/CPIC_Guideline_for_clopidogrel_and_CYP2C19.json"));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		Gson gson = new Gson();
		Rec recObj = gson.fromJson(br, Rec.class);
		System.out.println("ID Number: "+recObj.getId());
	}
}

