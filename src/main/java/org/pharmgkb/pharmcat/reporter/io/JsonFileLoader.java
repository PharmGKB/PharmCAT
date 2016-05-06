package org.pharmgkb.pharmcat.reporter.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import org.pharmgkb.pharmcat.haplotype.model.json.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.json.HaplotyperResult;
import org.pharmgkb.pharmcat.reporter.model.CPICException;
import org.pharmgkb.pharmcat.reporter.model.CPICExceptionList;
import org.pharmgkb.pharmcat.reporter.model.CPICinteraction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class JsonFileLoader {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  private final Gson gson = new Gson();

  public List<GeneCall> loadHaplotypeGeneCalls(@Nonnull File haplotypeCalledFile) throws IOException{
    Preconditions.checkNotNull(haplotypeCalledFile);
    Preconditions.checkArgument(haplotypeCalledFile.exists());
    Preconditions.checkArgument(haplotypeCalledFile.isFile());

    sf_logger.debug("Loading haplotyper file {}", haplotypeCalledFile);
    try (BufferedReader br = new BufferedReader(new FileReader( haplotypeCalledFile ))) {
      HaplotyperResult calls = gson.fromJson(br, HaplotyperResult.class);
      return calls.getGeneCalls();
    }
  }

  /**
   * This turns the json exception list into a fragmented list under the header
   * of the relevant gene.
   * FIXME dumb assumption that the exception are for single gene issues and not for
   * multi gene interactions
   */
  public Map<String, List<CPICException>> loadExceptions( File exceptions )throws IOException {
    Map<String, List<CPICException>> matcher = new HashMap<>();
    BufferedReader br = new BufferedReader(new FileReader( exceptions ));
    CPICExceptionList except = gson.fromJson(br, CPICExceptionList.class);
    for (CPICException rule : except.getRules()) {
      if( matcher.containsKey(rule.getGene())){
        matcher.get(rule.getGene()).add(rule);
      } else {
        List<CPICException> parts = new ArrayList<>();
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
  public List<CPICinteraction> loadDrugGeneRecommendations( List<File> interactions ) throws IOException {
    List<CPICinteraction> drugGenes = new ArrayList<>();
    for( File interact: interactions ){

      BufferedReader br = new BufferedReader(new FileReader( interact ));
      CPICinteraction act = gson.fromJson(br, CPICinteraction.class);
      drugGenes.add(act);
      br.close();
    }

    return drugGenes;
  }
}
