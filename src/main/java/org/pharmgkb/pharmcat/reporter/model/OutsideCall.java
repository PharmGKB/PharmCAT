package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.reporter.BadOutsideCallException;


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
  private static final int IDX_PHENO = 2;

  private String m_gene;
  private List<String> m_diplotypes;
  private String m_phenotype;
  private final SortedSet<String> m_haplotypes = new TreeSet<>(HaplotypeNameComparator.getComparator());

  /**
   * Constructor that expects a TSV formatted String with the gene symbol in the first field, the diplotype call
   * (e.g. *1/*2) in the second field, and the phenotype call in the third field.
   * @param line a TSV-formatted string
   * @throws RuntimeException can occur if data not in expected format
   */
  public OutsideCall(String line) throws RuntimeException {
    String[] fields = line.split(sf_lineSeparator);
    if (fields.length < 2) {
      throw new BadOutsideCallException("Expected at least 2 TSV fields, got " + fields.length);
    }

    String gene = fields[IDX_GENE];
    String dips = fields[IDX_DIPS];

    if (StringUtils.isBlank(gene)) {
      throw new BadOutsideCallException("No gene specified");
    }
    m_gene = gene;

    if (fields.length == 2 && StringUtils.isBlank(dips)) {
      throw new BadOutsideCallException("No diplotypes specified");
    }
    String diplotypes = dips.replaceAll(gene, "");
    m_diplotypes = new ArrayList<>();
    Arrays.stream(diplotypes.split(sf_separator))
        .map(StringUtils::trimToNull)
        .filter(Objects::nonNull)
        .forEach(d -> m_diplotypes.add(d));

    m_diplotypes.forEach(d -> {
      String[] alleles = d.split(sf_dipSeparator);
      if (alleles.length > 2) {
        throw new BadOutsideCallException("Too many alleles specified in " + d);
      }
      m_haplotypes.add(alleles[0]);
      if (alleles.length == 2) {
        m_haplotypes.add(alleles[1]);
      }
    });

    if (fields.length >= 3) {
      String pheno = fields[IDX_PHENO];
      m_phenotype = StringUtils.stripToNull(pheno.replaceAll(gene, ""));
      if (StringUtils.isBlank(m_phenotype) && m_diplotypes.size() == 0) {
        throw new BadOutsideCallException("Specify a diplotype or phenotype for " + gene);
      }
      if (StringUtils.isNotBlank(m_phenotype) && m_diplotypes.size() > 0) {
        throw new BadOutsideCallException("Specify either diplotype OR phenotype, not both for " + gene);
      }
    }
  }

  @Override
  public String toString() {
    if (getDiplotypes().size() > 0) {
      return getGene() + ":" + String.join(",", getDiplotypes());
    } else if (getPhenotype() != null) {
      return getGene() + ":" + getPhenotype();
    } else {
      return getGene();
    }
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

  public String getPhenotype() {
    return m_phenotype;
  }

  public void setPhenotype(String phenotype) {
    m_phenotype = phenotype;
  }
}
