package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Arrays;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.BinaryOperator;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.haplotype.model.Variant;


/**
 * Model class to represent a diplotype and all derived information
 *
 * @author Ryan Whaley
 */
public class Diplotype implements Comparable<Diplotype> {

  private static final String NA = "N/A";
  private static final String sf_toStringPattern = "%s:%s";
  private static final String sf_delimiter = "/";
  private static final String sf_homTemplate = "Two %s alleles";
  private static final String sf_hetTemplate = "One %s allele and one %s allele";

  private Haplotype m_allele1;
  private Haplotype m_allele2;
  private String m_gene;
  private String m_phenotype;
  private Variant m_variant;

  /**
   * This Function can be used in reduce() calls
   */
  public static final BinaryOperator<String> phasedReducer = (a,b)-> {
    if (a == null) {
      return b;
    }
    if (b == null) {
      return a;
    }

    int left = 0;
    int right = 1;

    String[] aHaps = a.split(sf_delimiter);
    String[] bHaps = b.split(sf_delimiter);

    Set<String> finalLeft = new TreeSet<>(HaplotypeNameComparator.getComparator());
    finalLeft.addAll(Arrays.asList(aHaps[left].split("\\+")));
    Set<String> finalRight = new TreeSet<>(HaplotypeNameComparator.getComparator());
    finalRight.addAll(Arrays.asList(aHaps[right].split("\\+")));

    if (finalLeft.contains(bHaps[right]) || finalRight.contains(bHaps[left])) {
      finalLeft.add(bHaps[right]);
      finalRight.add(bHaps[left]);
    }
    else {
      finalLeft.add(bHaps[left]);
      finalRight.add(bHaps[right]);
    }

    String joined0 = finalLeft.stream().collect(Collectors.joining("+"));
    String joined1 = finalRight.stream().collect(Collectors.joining("+"));

    return joined0+sf_delimiter+joined1;
  };

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

    String[] alleles = new String[]{m_allele1.printDisplay(), m_allele2.printDisplay()};
    Arrays.sort(alleles, HaplotypeNameComparator.getComparator());
    return Arrays.stream(alleles).collect(Collectors.joining(sf_delimiter));
  }

  /**
   * Gets a String representation of this Diplotype that can be used to display in output. This should <em>NOT</em> be
   * used for matching to annotation groups
   * @return a String display for this diplotype, without Gene symbol
   */
  public String printDisplay() {
    if (getVariant() != null) {
      return getVariant().printDisplay();
    }
    else {
      return printBare();
    }
  }

  /**
   * Gets a string diplotype pair used to look up matching guideline groups. This could be different than what's displayed in the
   * report to the user so we use a separate method.
   * @return a String key used to match guideline groups without gene symbol (e.g. *4/*10)
   */
  public String printBareLookupKey() {

    switch (m_gene) {
      case "DPYD":
        // this is here since *2B has no function assigned to it, remove once it's assigned
        if (getAllele1().getName().equals("*1") && getAllele2().getName().equals("*2B")) {
          return "*2A"+sf_delimiter+"*5";
        }
    }

    String[] alleles = new String[]{m_allele1.printLookup(), m_allele2.printLookup()};
    Arrays.sort(alleles, HaplotypeNameComparator.getComparator());
    return Arrays.stream(alleles).collect(Collectors.joining(sf_delimiter));
  }

  /**
   * Gets a string key used to look up matching guideline groups. This could be different than what's displayed in the
   * report to the user so we use a separate method.
   * @return a String key used to match guideline groups without gene symbol (e.g. *4/*10)
   */
  public String printLookupKey() {
    return m_gene + ":" + printBareLookupKey();
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
        if (getAllele1().isIncidental() && getAllele2().isIncidental()) {
          return Optional.of(getAllele1().getName() + " (homozygous)");
        }

        boolean override1 = getAllele1().getName().equals("Reference");
        boolean override2 = getAllele2().getName().equals("Reference");

        if (override1 && override2) {
          return Optional.of("No CPIC variants found");
        }
        else if (override1 || override2) {
          String allele = override1 ? getAllele2().getName() : getAllele1().getName();
          return Optional.of(allele + " (heterozygous)");
        }
        break;

    }
    return Optional.empty();
  }

  /**
   * Gets a variant used to make this diplotype call
   * @return a Variant used in this call
   */
  public Variant getVariant() {
    return m_variant;
  }

  public void setVariant(Variant variant) {
    m_variant = variant;
  }

  @Override
  public int compareTo(@Nonnull Diplotype o) {
    int rez = ObjectUtils.compare(getGene(), o.getGene());
    if (rez != 0) {
      return rez;
    }

    rez = ObjectUtils.compare(getAllele1(), o.getAllele1());
    if (rez != 0) {
      return rez;
    }

    return ObjectUtils.compare(getAllele2(), o.getAllele2());
  }
}
