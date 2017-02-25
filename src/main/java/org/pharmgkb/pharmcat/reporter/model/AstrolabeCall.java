package org.pharmgkb.pharmcat.reporter.model;

import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.annotation.Nonnull;
import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;


/**
 * Call data from astrolabe
 *
 * @author Ryan Whaley
 */
public class AstrolabeCall {

  private static final String sf_separator = " or ";
  private static final String sf_lineSeparator = "\t";
  private static final String sf_dipSeparator = "/";

  private static final int IDX_GENE = 0;
  private static final int IDX_DIPS = 1;

  private String m_gene;
  private List<String> m_diplotypes;
  private SortedSet<String> m_haplotypes = new TreeSet<>(HaplotypeNameComparator.getComparator());

  public AstrolabeCall(@Nonnull String line) throws RuntimeException {
    String[] fields = line.split(sf_lineSeparator);
    if (fields.length < 2) {
      throw new RuntimeException("Astrolabe output not in expected format");
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
      m_haplotypes.add(alleles[0]);
      m_haplotypes.add(alleles[1]);
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
