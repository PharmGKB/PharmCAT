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

  private static final Path VCF_HEADER_PATH = PathUtils.getPathToResource("org/pharmgkb/pharmcat/util/VcfHeader.txt");
  private static final String TEST_PATH  = "org/pharmgkb/pharmcat/haplotype/";

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
    return Files.lines(VCF_HEADER_PATH)
        .filter(l -> l.startsWith("#"))
        .collect(Collectors.joining("\n"));
  }
}
