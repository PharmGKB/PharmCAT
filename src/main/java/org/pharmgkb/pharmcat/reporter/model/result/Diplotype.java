package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Arrays;
import java.util.Optional;
import java.util.stream.Collectors;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;


/**
 * Model class to represent a diplotype and all derived information
 *
 * @author Ryan Whaley
 */
public class Diplotype {

  private static final String NA = "N/A";
  private static final String sf_toStringPattern = "%s:%s";
  private static final String sf_delimiter = "/";
  private static final String sf_homTemplate = "Two %s alleles";
  private static final String sf_hetTemplate = "One %s allele and one %s allele";

  private Haplotype m_allele1;
  private Haplotype m_allele2;
  private String m_gene;
  private String m_phenotype;

  /**
   * public constructor
   */
  public Diplotype(String gene, Haplotype h1, Haplotype h2) {

    m_allele1 = h1;
    m_allele2 = h2;
    m_gene = gene;
  }

  /**
   * Gets the gene this diplotype is for
   * @return a HGNC gene symbol
   */
  public String getGene() {
    return m_gene;
  }

  /**
   * Gets the first {@link Haplotype} listed in this diplotype
   */
  public Haplotype getAllele1() {
    return m_allele1;
  }

  /**
   * Gets the second {@link Haplotype} listed in this diplotype
   */
  public Haplotype getAllele2() {
    return m_allele2;
  }

  /**
   * Does this diplotype have an allele with the given name
   * @param alleleName an allele name, e.g. "*10"
   * @return true if this diplotype contains an allele with the given name
   */
  public boolean hasAllele(String alleleName) {
    return (m_allele1 != null && m_allele1.getName().equals(alleleName))
        || (m_allele2 != null && m_allele2.getName().equals(alleleName));
  }

  /**
   * Flag for whether this diplotype includes an incidental finding
   * @return true if this diplotype includes an incidental finding allele
   */
  public boolean isIncidental() {
    return m_allele1.isIncidental() || m_allele2.isIncidental();
  }

  /**
   * Gets a Sting representation of this haplotype with no gene prefix (e.g. *1/*10)
   */
  public String printBare() {

    Optional<String> override = printOverride();
    if (override.isPresent()) {
      return override.get();
    }

    String[] alleles = new String[]{m_allele1.getName(), m_allele2.getName()};
    Arrays.sort(alleles, HaplotypeNameComparator.getComparator());
    return Arrays.stream(alleles).collect(Collectors.joining(sf_delimiter));
  }

  /**
   * Gets a String phrase describing the individual haplotype functions, e.g. "Two low function alleles"
   *
   * Will print a default N/A String if no functions exist
   */
  public String printFunctionPhrase() {

    String f1 = getAllele1() != null && getAllele1().getFunction() != null ?
        getAllele1().getFunction().toLowerCase() : null;
    String f2 = getAllele2() != null && getAllele2().getFunction() != null ?
        getAllele2().getFunction().toLowerCase() : null;

    if (StringUtils.isNotBlank(f1) && StringUtils.isNotBlank(f2)) {

      if (f1.equals(f2)) {
        return String.format(sf_homTemplate, f1);
      }
      else {
        return String.format(sf_hetTemplate, f1, f2);
      }

    }
    return GeneReport.NA;
  }

  /**
   * Gets a String representation of this diplotype with the gene as prefix, e.g. GENEX:*1/*2
   */
  public String toString() {
    return String.format(sf_toStringPattern, m_gene, printBare());
  }

  /**
   * Gets a String term for the overall phenotype of this Diplotype
   *
   * Will print a default N/A String if no phenotype exists
   */
  public String getPhenotype() {
    return m_phenotype == null ? NA : m_phenotype;
  }

  public void setPhenotype(String phenotype) {
    m_phenotype = phenotype;
  }

  /**
   * Print the overriding diplotype string if it exists
   * @return Optional diplotype string to override whatever the actual string would be
   */
  private Optional<String> printOverride() {
    switch (m_gene) {

      case "CFTR":
        if (getAllele1().getName().equals("Reference") && getAllele1().getName().equals("Reference")) {
          return Optional.of("No CPIC variants found");
        }
        else if (getAllele1().getName().equals("Reference") || getAllele1().getName().equals("Reference")) {
          String allele = getAllele1().getName().equals("Reference") ? getAllele2().getName() : getAllele1().getName();
          return Optional.of(allele + " (heterozygous)");
        }
        break;

    }
    return Optional.empty();
  }
}
