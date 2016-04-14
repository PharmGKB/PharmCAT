package org.cpic.reporter.development;

import com.google.gson.Gson;
import org.cpic.reporter.model.CPICException;
import org.cpic.reporter.model.CPICExceptionList;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;



/**
 * Created by lester on 4/13/16.
 */
public class TestExceptionParser {
  public static void main(String[] args) throws Exception {
    BufferedReader br = null;
    try {
      br = new BufferedReader(new FileReader("resources/cpic_exceptions/exceptions.json"));
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    Gson gson = new Gson();

    CPICExceptionList cpicExceptions = gson.fromJson(br, CPICExceptionList.class);
    System.out.println(cpicExceptions.getRules()  );

    for (CPICException rule : cpicExceptions.getRules()) {
      System.out.println(rule.getGene());
      System.out.println("\t" + rule.getException_type());
      System.out.println("\t" + rule.getMatches());
      System.out.println("\t" + rule.getMessage());
    }


  }


}

