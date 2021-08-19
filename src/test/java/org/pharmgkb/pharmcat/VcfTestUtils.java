package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.io.StringWriter;
import java.nio.file.Files;
import org.pharmgkb.common.util.PathUtils;


/**
 * Helper methods for testing VCF files
 *
 * @author Ryan Whaley
 */
public class VcfTestUtils {

  private static final String TEST_PATH  = "org/pharmgkb/pharmcat/haplotype/";
  private static final String VCF_HEADER = "##fileformat=VCFv4.1\n" +
      "##fileDate=2015-08-04\n" +
      "##source=IlluminaPlatinumGenomes, version: hg38_2.0.1\n" +
      "##reference=hg38\n" +
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
      "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
      "##SAMPLE=<ID=PharmCAT,Description=\"Synthetic PharmCAT test\">\n" +
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n";

  public static String writeVcf(String[] filesToInclude) {
    StringWriter writer = new StringWriter();
    try {
      writer.write(VCF_HEADER);

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
}
