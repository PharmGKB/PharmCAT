package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Comparator;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.base.Charsets;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.util.DataManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class gathers information about intervals on the chromosomal sequence where PharmCAT needs data to make its
 * calls. This will be a list of 1 or more Strings in the format:
 *
 * chrXX:YYYYY-ZZZZZ
 *
 * Where XX is the chromosome name (e.g. 1, 21, X), YYYYY is the minimum position, and ZZZZZ is the maximum posiiton.
 *
 * There can be more than one interval for a chromosome and some chromosomes may not be listed at all.
 */
public class ChromosomalIntervals {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  private SortedSet<String> m_intervals = new TreeSet<>(Comparator.naturalOrder());

  /**
   * The CLI for writing the intervals to a given file path
   * @param args command-line arguments
   */
  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("o", "output-file", "output internal file", true, "file-path");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      ChromosomalIntervals intervalFile = new ChromosomalIntervals();
      intervalFile.write(cliHelper.getValidFile("o", false));
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  /**
   * Constructur. Will read through the definition files on instantiation and gather the interval data.
   * @throws IOException can occur from bad definition data
   */
  private ChromosomalIntervals() throws IOException {
    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);

    for (String geneSymbol : definitionReader.getGenes()) {
      sf_logger.info(geneSymbol);

      DefinitionFile defFile = definitionReader.getDefinitionFile(geneSymbol);

      long min = Arrays.stream(defFile.getVariants())
          .map(VariantLocus::getVcfPosition)
          .min(Long::compareTo)
          .orElseThrow(() -> new RuntimeException("No minimum for " + geneSymbol + " interval"));

      long max = Arrays.stream(defFile.getVariants())
          .map(VariantLocus::getVcfPosition)
          .max(Long::compareTo)
          .orElseThrow(() -> new RuntimeException("No maximum for " + geneSymbol + " interval"));

      sf_logger.info(String.format("%s: %d-%d", defFile.getChromosome(), min, max));
      m_intervals.add(String.format("%s: %d-%d", defFile.getChromosome(), min, max));
    }
  }

  /**
   * Write the already gathered intervals to the given file
   */
  private void write(Path outputPath) throws IOException {
    Files.write(outputPath, m_intervals, Charsets.UTF_8);
    sf_logger.info("Wrote data to {}", outputPath);
  }
}
