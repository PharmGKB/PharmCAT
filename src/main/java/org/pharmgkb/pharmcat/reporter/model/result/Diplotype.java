package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import javax.annotation.Nonnull;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;


/**
 * Model class to represent a diplotype and all derived information
 *
 * @author Ryan Whaley
 */
public class Diplotype implements Comparable<Diplotype> {
  public static final String UNKNOWN = "Unknown/Unknown";

  private static final String NA = "N/A";
  private static final String sf_toStringPattern = "%s:%s";
  private static final String sf_delimiter = "/";
  private static final String sf_homTemplate = "Two %s alleles";
  private static final String sf_hetTemplate = "One %s allele and one %s allele";
  private static final String sf_hetSuffix = " (heterozygous)";
  private static final String sf_homSuffix = " (homozygous)";

  @Expose
  @SerializedName("allele1")
  private final Haplotype m_allele1;
  @Expose
  @SerializedName("allele2")
  private final Haplotype m_allele2;
  @Expose
  @SerializedName("gene")
  private final String m_gene;
  @Expose
  @SerializedName("phenotype")
  private String m_phenotype;
  @Expose
  @SerializedName("variant")
  private VariantReport m_variant;
  @Expose
  @SerializedName("lookupKey")
  private String m_lookupKey;
  @Expose
  @SerializedName("label")
  private final String f_label;

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
    Set<String> finalRight = new TreeSet<>(HaplotypeNameComparator.getComparator());

    addHaps(finalLeft, aHaps[left]);
    addHaps(finalRight, aHaps[right]);

    if (finalLeft.contains(bHaps[right]) || finalRight.contains(bHaps[left])) {
      addHaps(finalLeft, bHaps[right]);
      addHaps(finalRight, bHaps[left]);
    }
    else {
      addHaps(finalLeft, bHaps[left]);
      addHaps(finalRight, bHaps[right]);
    }

    String joined0 = String.join("+", finalLeft);
    String joined1 = String.join("+", finalRight);

    return joined0+sf_delimiter+joined1;
  };

  /**
   * When multiple diplotype calls are phased they can be combined into one String in certain circumstances. This method
   * combines things like a/b and a/c into a/b+c.
   * @param dips a {@link Collection} of diplotype strings in the form "a/b", "a/c", etc... 
   * @return a single String combining the diplotypes, e.g. "a/b+c"
   */
  @Nonnull
  public static String reducePhasedDiplotypes(@Nullable Collection<String> dips) {
    if (dips == null || dips.size() == 0) return "";
    
    final Set<String> leftBucket = new TreeSet<>(HaplotypeNameComparator.getComparator());
    final Set<String> rightBucket = new TreeSet<>(HaplotypeNameComparator.getComparator());
    
    dips.forEach(d -> {
      String[] haps = d.split(sf_delimiter);
      if (!rightBucket.contains(haps[0])) {
        leftBucket.add(haps[0]);
        rightBucket.add(haps[1]);
      } else {
        leftBucket.add(haps[1]);
        rightBucket.add(haps[0]);
      }
    });
    
    Function<String,String> leftMapper = leftBucket.size() > 1 ? hapNameEncloser : Function.identity();
    String left = leftBucket.stream().map(leftMapper).collect(Collectors.joining("+"));
    Function<String,String> rightMapper = rightBucket.size() > 1 ? hapNameEncloser : Function.identity();
    String right = rightBucket.stream().map(rightMapper).collect(Collectors.joining("+"));
    
    return left + sf_delimiter + right;
  }

  /**
   * Takes diplotype names and encloses them in parentheses if they contain a "+"
   */
  private static final Function<String,String> hapNameEncloser = h -> {
    if (h.contains("+")) return "(" + h + ")"; else return h;
  };

  private static void addHaps(Set<String> hapSet, String hap) {
    hapSet.addAll(Arrays.asList(hap.split("\\+")));
  }

  /**
   * public constructor
   */
  public Diplotype(String gene, Haplotype h1, Haplotype h2) {
    m_allele1 = h1;
    m_allele2 = h2;
    m_gene = gene;
    f_label = printBare();
  }

  public Diplotype(String gene, Haplotype h) {
    m_allele1 = h;
    m_allele2 = null;
    m_gene = gene;
    f_label = printBare();
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
  private Haplotype getAllele1() {
    return m_allele1;
  }

  /**
   * Gets the second {@link Haplotype} listed in this diplotype
   */
  private Haplotype getAllele2() {
    return m_allele2;
  }

  /**
   * Is this Diplotype for a single-ploidy gene like MT-RNR1. We acknowledge this isn't a "real" diplotype if this is
   * true but for the purposes of this system we will call it a diplotype.
   *
   * @return true if this Diplotype is single ploidy
   */
  public boolean isSinglePloidy() {
    return m_allele2 == null;
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

  public boolean isUnknown() {
    return m_allele1.isUnknown() && m_allele2.isUnknown();
  }

  /**
   * Gets a Sting representation of this haplotype with no gene prefix (e.g. *1/*10)
   */
  public String printBare() {
    return printOverride().orElseGet(() -> {
      if (m_allele2 != null) {
        String[] alleles = new String[]{ m_allele1.getName(), m_allele2.getName() };
        Arrays.sort(alleles, HaplotypeNameComparator.getComparator());
        return String.join(sf_delimiter, alleles);
      } else {
        return m_allele1.getName();
      }
    });
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
   * Gets a String phrase describing the individual haplotype functions, e.g. "Two low function alleles"
   *
   * Will print a default N/A String if no functions exist
   */
  public String printFunctionPhrase() {

    String f1 = getAllele1() != null && getAllele1().getFunction() != null ?
        getAllele1().getFunction().toLowerCase() : null;
    String f2 = getAllele2() != null && getAllele2().getFunction() != null ?
        getAllele2().getFunction().toLowerCase() : null;

    if (!isSinglePloidy() && StringUtils.isNotBlank(f1) && StringUtils.isNotBlank(f2)) {

      if (f1.equals(f2)) {
        return String.format(sf_homTemplate, f1);
      }
      else {
        String[] functions = new String[]{f1, f2};
        Arrays.sort(functions);
        return String.format(sf_hetTemplate, functions[0], functions[1]);
      }

    } else if (isSinglePloidy() && StringUtils.isNotBlank(f1)) {
      return f1;
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
   *
   * <strong>Use only for display purposes, not for matchign recommendations</strong>
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
    boolean refAllele1 = getAllele1().isReference();
    boolean refAllele2 = !isSinglePloidy() && getAllele2().isReference();

    if (m_gene.equals("CFTR")) {
      if (refAllele1 && refAllele2) {
        return Optional.of("No CPIC variants found");
      }
      else if (refAllele1 || refAllele2) {
        String allele = refAllele1 ? getAllele2().getName() : getAllele1().getName();
        return Optional.of(allele + sf_hetSuffix);
      }
    }
    return Optional.empty();
  }

  /**
   * Makes a {@link Stream} of zygosity descriptors for this diplotype, e.g. *60 (heterozygous). This is a stream since
   * a single Diplotype can be described in 0, 1, or 2 Strings, depending on the particular allele.
   * @return a Stream of 0 or more zygosity Strings
   */
  public Stream<String> streamAllelesByZygosity() {
    if (getAllele1().equals(getAllele2())) {
      if (getAllele1().isReference()) {
        return Stream.empty();
      }
      return Stream.of(getAllele1().getName() + sf_homSuffix);
    }
    else {
      Set<String> alleles = new TreeSet<>(HaplotypeNameComparator.getComparator());
      if (!getAllele1().isReference()) {
        alleles.add(getAllele1().getName() + sf_hetSuffix);
      }
      if (!getAllele2().isReference()) {
        alleles.add(getAllele2().getName() + sf_hetSuffix);
      }
      return alleles.stream();
    }
  }

  /**
   * Gets a variant used to make this diplotype call
   * @return a Variant used in this call
   */
  public VariantReport getVariant() {
    return m_variant;
  }

  public void setVariant(VariantReport variant) {
    m_variant = variant;
  }

  @Override
  public int compareTo(Diplotype o) {
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

  public String getLookupKey() {
    return m_lookupKey;
  }

  public void setLookupKey(String lookupKey) {
    m_lookupKey = lookupKey;
  }

  public Map<String,Integer> makeLookupMap() {
    Map<String,Integer> lookupMap = new HashMap<>();
    if (m_allele2 != null) {
      if (m_allele1.equals(m_allele2)) {
        lookupMap.put(m_allele1.getName(), 2);
      } else {
        lookupMap.put(m_allele1.getName(), 1);
        lookupMap.put(m_allele2.getName(), 1);
      }
    } else {
      lookupMap.put(m_allele1.getName(), 1);
    }
    return lookupMap;
  }
}
