package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

/**
 * This will give a map of gene and haplotype name to unique ID.
 *
 * @author Alex Frase
 * @author Ryan Whaley
 */
public class HaplotypeIdMap {
  private Map<String,Map<String,String>> m_geneHaplotypeMap = new HashMap<>();


  public HaplotypeIdMap() {
    try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("haplotype.id.list.tsv")))) {
      // first line is the sequence
      reader.readLine();
      String line;
      while ((line = reader.readLine()) != null) {
        String[] fields = line.split("\t");

        Map<String, String> hapMap = m_geneHaplotypeMap.computeIfAbsent(fields[0], k -> new HashMap<>());
        hapMap.put(fields[1], fields[2]);
      }
    } catch (IOException ex) {
      throw new IllegalStateException("Cannot find haplotype.id.list.tsv", ex);
    }
  }


  /**
   * Gets the unique ID for a given gene and haplotype identifier.
   * @param geneName a gene name (i.e. "CYP2C9")
   * @param haplotypeName a haplotype name (i.e. "*3")
   * @return a String representation of a unique ID
   */
  public String getId(String geneName, String haplotypeName) {
    if (!m_geneHaplotypeMap.containsKey(geneName)) {
      return null;
    }
    return m_geneHaplotypeMap.get(geneName).get(haplotypeName);
  }
}
