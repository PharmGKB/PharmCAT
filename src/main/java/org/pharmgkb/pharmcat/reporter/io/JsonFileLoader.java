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
import org.pharmgkb.pharmcat.haplotype.Haplotyper;
import org.pharmgkb.pharmcat.haplotype.model.json.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.json.HaplotyperResult;
import org.pharmgkb.pharmcat.reporter.model.CPICinteraction;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class takes JSON files and deserializes them into Objects through GSON.
 */
public class JsonFileLoader {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  private final Gson gson = new Gson();

  /**
   * Load all the gene calls coming from the {@link Haplotyper} utility
   */
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
   * Load all the guideline annotations into {@link CPICinteraction} objects from the list of guideline {@link File}
   * list
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
