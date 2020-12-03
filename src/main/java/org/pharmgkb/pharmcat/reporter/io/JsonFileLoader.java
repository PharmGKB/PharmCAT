package org.pharmgkb.pharmcat.reporter.io;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import com.google.common.base.Preconditions;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class takes JSON files and deserializes them into Objects through GSON.
 */
public class JsonFileLoader {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

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
}
