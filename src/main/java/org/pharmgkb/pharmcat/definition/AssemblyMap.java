package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;
import org.pharmgkb.common.util.PathUtils;


/**
 * This will give a map of RefSeq chromosome identifier to GRC build number.
 *
 * @author Ryan Whaley
 */
public class AssemblyMap {
  public static final String GRCH38 = "b38";
  public static final String GRCH37 = "b37";
  private static final Path sf_mappingFile =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/chr_build_mapping.tsv");
  private final Map<String,String> m_assemblyMap = new HashMap<>();

  public AssemblyMap() throws IOException {
    try (BufferedReader reader = Files.newBufferedReader(sf_mappingFile)) {
      String line;
      while ((line = reader.readLine()) != null) {
        String[] fields = line.split("\t");
        m_assemblyMap.put(fields[3], fields[0]);
      }
    }
  }

  /**
   * Gets the build number (in format b##) for a given RefSeq chromosome identifier (e.g. NC_00001.10)
   * @param refSeqId a RefSeq chromosome identifier
   * @return a String representation of a build number
   */
  public String get(String refSeqId) {
    return m_assemblyMap.get(refSeqId);
  }
}
