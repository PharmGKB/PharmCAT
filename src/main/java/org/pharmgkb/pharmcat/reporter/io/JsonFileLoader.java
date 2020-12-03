package org.pharmgkb.pharmcat.reporter.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.reporter.model.DosingGuideline;
import org.pharmgkb.pharmcat.reporter.model.GuidelinePackage;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class takes JSON files and deserializes them into Objects through GSON.
 */
public class JsonFileLoader {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final String CPIC_SOURCE = "CPIC";

  private final Gson gson = new Gson();

  /**
   * Load all the gene calls coming from the {@link NamedAlleleMatcher} utility
   */
  public List<GeneCall> loadHaplotypeGeneCalls(Path haplotypeCalledFile) throws IOException{
    Preconditions.checkNotNull(haplotypeCalledFile);
    Preconditions.checkArgument(Files.exists(haplotypeCalledFile));
    Preconditions.checkArgument(Files.isRegularFile(haplotypeCalledFile));

    sf_logger.debug("Loading haplotyper file {}", haplotypeCalledFile);
    Result namResult = new ResultSerializer().fromJson(haplotypeCalledFile);
    return namResult.getGeneCalls();
  }

  /**
   * Load the <strong>CPIC</strong> guideline annotations into {@link DosingGuideline} objects from the list of guideline {@link File}
   * list
   * @deprecated replace with drug list
   */
  public List<GuidelinePackage> loadGuidelines(List<Path> guidelineFileList) throws IOException {
    List<GuidelinePackage> guidelines = new ArrayList<>();

    for (Path guidelineFile : guidelineFileList) {
      try (BufferedReader br = Files.newBufferedReader(guidelineFile)) {
        GuidelinePackage guidelinePackage = gson.fromJson(br, GuidelinePackage.class);
        DosingGuideline guideline = guidelinePackage.getGuideline();

        if (guideline.getSource().equals(CPIC_SOURCE) && guideline.isRecommendation()) {
          guidelines.add(guidelinePackage);
        }
      }
    }

    return guidelines;
  }
}
