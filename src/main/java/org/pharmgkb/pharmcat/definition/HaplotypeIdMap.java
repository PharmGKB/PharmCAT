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
  private int m_idSequence;
  private Map<String,Map<String,String>> m_geneHaplotypeMap = new HashMap<>();


  public HaplotypeIdMap() throws IOException {
    try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("haplotype.id.list.tsv")))) {
      // first line is the sequence
      m_idSequence = Integer.parseInt(reader.readLine());
      String line;
      while ((line = reader.readLine()) != null) {
        String[] fields = line.split("\t");

        Map<String, String> hapMap = m_geneHaplotypeMap.get(fields[0]);
        if (hapMap == null) {
          hapMap = new HashMap<>();
          m_geneHaplotypeMap.put(fields[0], hapMap);
        }
        hapMap.put(fields[1], fields[2]);
      }
    }
  }


  public String newId(String geneName, String haplotypeName) {

    if (m_geneHaplotypeMap.containsKey(geneName) && m_geneHaplotypeMap.get(geneName).containsKey(haplotypeName)) {
      throw new IllegalStateException(geneName + " " + haplotypeName + " has already been assigned ID "+
          m_geneHaplotypeMap.get(geneName).get(haplotypeName));
    }

    String id = (m_idSequence + 1) + ".1";
    Map<String, String> hapMap = m_geneHaplotypeMap.get(geneName);
    if (hapMap == null) {
      hapMap = new HashMap<>();
      m_geneHaplotypeMap.put(geneName, hapMap);
    }
    hapMap.put(haplotypeName, id);
    m_idSequence += 1;
    return id;
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
