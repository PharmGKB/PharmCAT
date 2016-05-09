package org.pharmgkb.pharmcat.reporter.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.List;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import org.pharmgkb.pharmcat.haplotype.model.json.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.json.HaplotyperResult;
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
   * DANGER WILL ROBINSON! A DUMB ASSUMPTION HAS BEEN MADE!
   * FIXME : assumption of single gene and single drug interaction here
   *
   * Multi gene drugs need alteration in this file to work but this will be
   * acceptable for hackathon standards but should be fixed as part of a the
   * next version.  This is a critical flaw in the design
   *
   */
  public List<CPICinteraction> loadGuidelines(List<File> guidelineFileList) throws IOException {
    List<CPICinteraction> drugGenes = new ArrayList<>();

    for (File guideline : guidelineFileList) {
      try (BufferedReader br = new BufferedReader(new FileReader(guideline))) {
        drugGenes.add(gson.fromJson(br, CPICinteraction.class));
      }
    }

    return drugGenes;
  }
}
