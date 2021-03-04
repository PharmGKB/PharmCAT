package org.pharmgkb.pharmcat.reporter.model;

import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.ParseException;


/**
 * Calls from an outside data source. For example, Astrolabe (but could be others)
 *
 * @author Ryan Whaley
 */
public class OutsideCall {

  private static final String sf_separator = " or ";
  private static final String sf_lineSeparator = "\t";
  private static final String sf_dipSeparator = "/";

  private static final int IDX_GENE = 0;
  private static final int IDX_DIPS = 1;

  private String m_gene;
  private List<String> m_diplotypes;
  private final SortedSet<String> m_haplotypes = new TreeSet<>(HaplotypeNameComparator.getComparator());

  /**
   * Constructor that expects a TSV formatted String with the gene symbol in the first field and the diplotype call 
   * (e.g. *1/*2) in the second field.
   * @param line a TSV-formatted string
   * @throws RuntimeException can occur if data not in expected format
   */
  public OutsideCall(String line) throws RuntimeException {
    String[] fields = line.split(sf_lineSeparator);
    if (fields.length < 2) {
      throw new RuntimeException("Expected at least 2 TSV fields, got " + fields.length);
    }

    String gene = fields[IDX_GENE];
    if (StringUtils.isBlank(gene)) {
      throw new RuntimeException("No gene specified");
    }
    m_gene = gene;

    String dips = fields[IDX_DIPS];
    if (StringUtils.isBlank(dips)) {
      throw new RuntimeException("No diplotypes specified");
    }
    String diplotypes = dips.replaceAll(gene, "");
    m_diplotypes = ImmutableList.copyOf(diplotypes.split(sf_separator));

    m_diplotypes.forEach(d -> {
      String[] alleles = d.split(sf_dipSeparator);
      if (alleles.length > 2) {
        throw new ParseException("Too many alleles specified in " + d);
      }
      m_haplotypes.add(alleles[0]);
      if (alleles.length == 2) {
        m_haplotypes.add(alleles[1]);
      }
    });
  }

  @Override
  public String toString() {
    return getGene();
  }

  public String getGene() {
    return m_gene;
  }

  public void setGene(String gene) {
    m_gene = gene;
  }

  public List<String> getDiplotypes() {
    return m_diplotypes;
  }

  public void setDiplotypes(List<String> diplotypes) {
    m_diplotypes = diplotypes;
  }

  public SortedSet<String> getHaplotypes() {
    return m_haplotypes;
  }
}
