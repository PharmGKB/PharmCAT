package org.pharmgkb.pharmcat.phenotype.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.NonNull;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.ComparisonChain;
import org.pharmgkb.pharmcat.Constants;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.haplotype.CombinationMatcher;
import org.pharmgkb.pharmcat.phenotype.PhenotypeUtils;
import org.pharmgkb.pharmcat.reporter.BadOutsideCallException;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * Calls from an outside data source.
 *
 * @author Ryan Whaley
 */
public class OutsideCall implements Comparable<OutsideCall> {
  public static final Pattern CYP2D6_SUBALLELE_PATTERN = Pattern.compile("(\\*\\d+)\\.\\d+");
  public static final Pattern HLA_SUBALLELE_PATTERN = Pattern.compile("(\\*\\d+:\\d+)((?::\\d+)+)");
  private static final Splitter sf_lineSplitter = Splitter.on("\t").trimResults();
  private static final Splitter sf_diplotypeSplitter = Splitter.on("/").trimResults();
  private static final String sf_diplotypeSeparator = "/";

  private static final int IDX_GENE = 0;
  private static final int IDX_DIPS = 1;
  private static final int IDX_PHENO = 2;
  private static final int IDX_ACTIVITY = 3;

  private final String m_gene;
  private @Nullable String m_diplotype;
  private List<String> m_diplotypes;
  private @Nullable String m_phenotype = null;
  private @Nullable String m_activityScore = null;
  private final SortedSet<String> m_haplotypes = new TreeSet<>(HaplotypeNameComparator.getComparator());
  private final List<String> m_warnings = new ArrayList<>();

  /**
   * Constructor that expects a TSV formatted String with the gene symbol in the first field, the diplotype call
   * (e.g. *1/*2) in the second field, and the phenotype call in the third field.
   * @param line a TSV-formatted string
   * @throws RuntimeException can occur if data not in expected format
   */
  public OutsideCall(Env env, String line, int lineNumber) throws RuntimeException {
    List<String> fields = sf_lineSplitter.splitToList(line);
    if (fields.size() < 2) {
      throw new BadOutsideCallException("Line " + lineNumber + ": Expected at least 2 TSV fields, got " + fields.size());
    }

    m_gene = StringUtils.stripToNull(fields.get(IDX_GENE));
    if (m_gene == null) {
      throw new BadOutsideCallException("Line " + lineNumber + ": No gene specified");
    }
    m_diplotype = StringUtils.stripToNull(fields.get(IDX_DIPS));
    if (fields.size() >= 3) {
      m_phenotype = PhenotypeUtils.normalize(fields.get(IDX_PHENO).replaceAll(m_gene, ""));
    }
    if (fields.size() >= 4) {
      m_activityScore = StringUtils.stripToNull(fields.get(IDX_ACTIVITY));
    }

    if (!env.hasGene(m_gene)) {
      return;
    }
    if (m_phenotype == null && m_activityScore == null && m_diplotype == null) {
      throw new BadOutsideCallException("Specify a diplotype, phenotype, or activity score for " + m_gene);
    }

    if (fields.size() == 2 && (m_diplotype == null || m_diplotype.equals(sf_diplotypeSeparator))) {
      if (StringUtils.isBlank(fields.get(IDX_DIPS))) {
        throw new BadOutsideCallException("Line " + lineNumber + ": No diplotype specified");
      } else {
        throw new BadOutsideCallException("Line " + lineNumber + ": Invalid diplotype specified");
      }
    }

    if (m_diplotype == null) {
      m_diplotypes = ImmutableList.of();
      if (env.isActivityScoreGene(m_gene)) {
        m_warnings.add(m_gene + " is not an activity score gene but has outside call with only an " +
            "activity score.  PharmCAT will not be able to provide any recommendations based on this gene.");
      }
    } else {
      // strip any prefix of the gene symbol
      List<String> alleles = sf_diplotypeSplitter.splitToList(m_diplotype).stream()
          .map(a -> a.replaceFirst("^" + m_gene + "\\s*", ""))
          .toList();
      if (alleles.size() > 2) {
        throw new BadOutsideCallException("Line " + lineNumber + ": Too many alleles specified in " + m_diplotype);
      }

      if (m_gene.equals("CYP2D6")) {
        alleles = alleles.stream()
            .map(a -> {
              Matcher m = CYP2D6_SUBALLELE_PATTERN.matcher(a);
              if (m.matches()) {
                m_warnings.add("PharmCAT does not support sub-alleles for " + m_gene + ". Using '" + m.group(1) +
                    "' instead of '" + a + "'.");
                return m.group(1);
              }
              return a;
            })
            .toList();

      } else if (m_gene.equals("HLA-A") || m_gene.equals("HLA-B")) {
        String prefix = m_gene.substring(m_gene.length() - 1) + "*";
        alleles = alleles.stream()
            .map(a -> {
              String orig = a;
              if (!a.startsWith("*")) {
                if (!a.startsWith(prefix)) {
                  throw new BadOutsideCallException("Invalid " + m_gene + " allele: '" + orig + "'.");
                }
                a = a.substring(1);
              }
              Matcher m = HLA_SUBALLELE_PATTERN.matcher(a);
              boolean isSuballele = false;
              if (m.matches()) {
                a = m.group(1);
                isSuballele = true;
              }
              if (!orig.equals(a)) {
                if (isSuballele) {
                  m_warnings.add("PharmCAT does not support sub-alleles for " + m_gene + ". Using '" + a +
                      "' instead of '" + orig + "'.");
                } else {
                  m_warnings.add("Converting outside call for " + m_gene + " from '" + orig + "', to '" + a + "'.");
                }
              }
              return a;
            })
            .toList();
      }

      // convert alleles from combination format if applicable
      alleles = alleles.stream()
          .map(a -> {
            if (!env.isValidNamedAllele(m_gene, a)) {
              String fixedA;
              if (a.startsWith("[") && a.endsWith("]")) {
                // convert PharmCAT style combinations into combinations recognized by phenotyper
                fixedA = a.substring(1, a.length() - 1);
                if (env.isValidNamedAllele(m_gene, fixedA)) {
                  m_warnings.add("Converting outside call for " + m_gene + " from '" + a + "', to '" + fixedA +
                      "'.");
                  return fixedA;
                }
              } else {
                fixedA = a;
              }
              if (fixedA.contains(CombinationMatcher.COMBINATION_JOINER)) {
                fixedA = fixedA.replaceAll(CombinationMatcher.COMBINATION_JOINER_REGEX, "+");
                if (env.isValidNamedAllele(m_gene, fixedA)) {
                  m_warnings.add("Converting outside call for " + m_gene + " from '" + a + "', to '" + fixedA +
                      "'.");
                  return fixedA;
                }
              }
              StringBuilder builder = new StringBuilder().append("Undocumented ")
                  .append(m_gene)
                  .append(" named ");
              if (Constants.isVariantGene(m_gene)) {
                builder.append("variant");
              } else {
                builder.append("allele");
              }
              builder.append(" in outside call: ")
                  .append(a);
              m_warnings.add(builder.toString());
            }
            return a;
          })
          .toList();

      // re-join alleles to eliminate white space when a gene symbol is used in diplotype
      m_diplotype = String.join(sf_diplotypeSeparator, alleles);
      m_diplotypes = ImmutableList.of(m_diplotype);

      m_haplotypes.add(alleles.get(0));
      if (alleles.size() == 2) {
        m_haplotypes.add(alleles.get(1));
      }
    }

    if (fields.size() >= 3) {
      m_phenotype = PhenotypeUtils.normalize(fields.get(IDX_PHENO).replaceAll(m_gene, ""));
    }

    if (fields.size() >= 4) {
      m_activityScore = StringUtils.stripToNull(fields.get(IDX_ACTIVITY));
    }
  }

  public void addWarning(String warning) {
    m_warnings.add(warning);
  }

  @Override
  public String toString() {
    if (m_diplotype != null) {
      return m_gene + ":" + m_diplotype;
    } else if (m_phenotype != null) {
      return m_gene + ":" + m_phenotype;
    }
    return m_gene + ":" + m_activityScore;
  }

  public String getGene() {
    return m_gene;
  }

  public @Nullable String getDiplotype() {
    return m_diplotype;
  }

  public List<String> getDiplotypes() {
    return m_diplotypes;
  }

  public SortedSet<String> getHaplotypes() {
    return m_haplotypes;
  }

  public @Nullable String getPhenotype() {
    return m_phenotype;
  }

  public @Nullable String getActivityScore() {
    return m_activityScore;
  }


  public List<String> getWarnings() {
    return m_warnings;
  }


  @Override
  public boolean equals(Object o) {
    if (!(o instanceof OutsideCall oc)) {
      return false;
    }
    if (o == this) {
      return true;
    }
    // ignore diplotypes and haplotypes because they are derived props
    return Objects.equals(m_gene, oc.getGene()) &&
        Objects.equals(m_diplotype, oc.getDiplotype()) &&
        Objects.equals(m_phenotype, oc.getPhenotype()) &&
        Objects.equals(m_activityScore, oc.getActivityScore())
        ;
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_gene, m_diplotype, m_phenotype, m_activityScore);
  }

  @Override
  public int compareTo(@NonNull OutsideCall o) {
    if (this == o) {
      return 0;
    }
    // ignore diplotypes and haplotypes because they are derived props
    return new ComparisonChain()
        .compare(m_gene, o.getGene())
        .compare(m_diplotype, o.getDiplotype())
        .compare(m_phenotype, o.getPhenotype())
        .compare(m_activityScore, o.getActivityScore())
        .result();
  }
}
