package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.stream.Collectors;
import org.pharmgkb.common.util.PathUtils;


/**
 * Helper methods for testing VCF files
 *
 * @author Ryan Whaley
 */
public class VcfTestUtils {

  private static final String HEADER_FILE = "org/pharmgkb/pharmcat/haplotype/DPYD/s1s1.vcf";
  private static final String TEST_PATH   = "org/pharmgkb/pharmcat/haplotype/";

  public static String writeVcf(String[] filesToInclude) {
    StringWriter writer = new StringWriter();
    try {
      writer.write(getVcfHeader());
      writer.write("\n");

      for (String filepath : filesToInclude) {
        Files.lines(PathUtils.getPathToResource(TEST_PATH + filepath))
            .filter(l -> !l.startsWith("#"))
            .forEach(l -> writer.write(l + "\n"));
      }
    } catch (IOException e) {
      e.printStackTrace();
    }

    return writer.toString();
  }

  private static String getVcfHeader() throws IOException {
    Path headerFile = PathUtils.getPathToResource(HEADER_FILE);

    return Files.lines(headerFile)
        .filter(l -> l.startsWith("#"))
        .collect(Collectors.joining("\n"));
  }
}
