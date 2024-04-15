package org.pharmgkb.pharmcat.phenotype.model;

import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.NonNull;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.ComparisonChain;
import org.pharmgkb.pharmcat.phenotype.PhenotypeUtils;
import org.pharmgkb.pharmcat.reporter.BadOutsideCallException;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * Calls from an outside data source.
 *
 * @author Ryan Whaley
 */
public class OutsideCall implements Comparable<OutsideCall> {
  private static final Splitter sf_lineSplitter = Splitter.on("\t").trimResults();
  private static final Splitter sf_diplotypeSplitter = Splitter.on("/").trimResults();
  private static final String sf_diplotypeSeparator = "/";

  private static final int IDX_GENE = 0;
  private static final int IDX_DIPS = 1;
  private static final int IDX_PHENO = 2;
  private static final int IDX_ACTIVITY = 3;

  private final String m_gene;
  private final @Nullable String m_diplotype;
  private final List<String> m_diplotypes;
  private @Nullable String m_phenotype = null;
  private @Nullable String m_activityScore = null;
  private final SortedSet<String> m_haplotypes = new TreeSet<>(HaplotypeNameComparator.getComparator());

  /**
   * Constructor that expects a TSV formatted String with the gene symbol in the first field, the diplotype call
   * (e.g. *1/*2) in the second field, and the phenotype call in the third field.
   * @param line a TSV-formatted string
   * @throws RuntimeException can occur if data not in expected format
   */
  public OutsideCall(String line, int lineNumber) throws RuntimeException {
    List<String> fields = sf_lineSplitter.splitToList(line);
    if (fields.size() < 2) {
      throw new BadOutsideCallException("Line " + lineNumber + ": Expected at least 2 TSV fields, got " + fields.size());
    }

    m_gene = StringUtils.stripToNull(fields.get(IDX_GENE));
    if (m_gene == null) {
      throw new BadOutsideCallException("Line " + lineNumber + ": No gene specified");
    }

    String diplotype = StringUtils.stripToNull(fields.get(IDX_DIPS));
    if (fields.size() == 2 && (diplotype == null || diplotype.equals(sf_diplotypeSeparator))) {
      if (StringUtils.isBlank(fields.get(IDX_DIPS))) {
        throw new BadOutsideCallException("Line " + lineNumber + ": No diplotype specified");
      } else {
        throw new BadOutsideCallException("Line " + lineNumber + ": Invalid diplotype specified");
      }
    }

    if (diplotype == null) {
      m_diplotype = null;
      m_diplotypes = ImmutableList.of();
    } else {
      // strip any prefix of the gene symbol
      List<String> alleles = sf_diplotypeSplitter.splitToList(diplotype).stream()
          .map(a -> a.replaceFirst("^" + m_gene + "\\s*", ""))
          .toList();
      // re-join alleles to eliminate white space when gene symbol is used in diplotype
      m_diplotype = String.join(sf_diplotypeSeparator, alleles);
      m_diplotypes = ImmutableList.of(m_diplotype);
      if (alleles.size() > 2) {
        throw new BadOutsideCallException("Line " + lineNumber + ": Too many alleles specified in " + m_diplotype);
      }
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

    if (m_phenotype == null && m_activityScore == null && m_diplotype == null) {
      throw new BadOutsideCallException("Specify a diplotype, phenotype, or activity score for " + m_gene);
    }
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
