package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import com.google.gson.Gson;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;


/**
 * This class is to help detect "incidental finding" alleles. It reads in the current list of incidental alleles from
 * and included file and exposes a method to detect whether a particular allele in an incidental finding.
 *
 * @author Ryan Whaley
 */
public class IncidentalFinder {

  private static final String INCIDENTAL_FILE = "incidental.alleles.json";

  private Map m_geneAlleleMap;

  /**
   * public constructor
   *
   * initializes incidental finding list from a file in the codebase
   */
  public IncidentalFinder() throws Exception {
    try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(INCIDENTAL_FILE)))) {

      Gson gson = new Gson();
      m_geneAlleleMap = gson.fromJson(reader, Map.class);

    } catch (IOException e) {
      throw new Exception("Error reading phenotype definitions", e);
    }
  }

  /**
   * Checks whether the given haplotype is an incidental finding
   *
   * @param haplotype a {@link Haplotype} object
   * @return true if the haplotype name is a known incidental finding
   */
  public boolean isFinding(Haplotype haplotype) {
    if (haplotype == null) {
      return false;
    }

    List alleles = (List)m_geneAlleleMap.get(haplotype.getGene());

    return alleles != null && alleles.contains(haplotype.getName());
  }

}
