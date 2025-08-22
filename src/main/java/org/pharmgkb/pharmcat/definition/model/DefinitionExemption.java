package org.pharmgkb.pharmcat.definition.model;

import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;


/**
 * This class represents special exemptions that must be applied to the {@link DefinitionFile} for the given
 * {@code gene}.
 *
 * @author Mark Woon
 */
public class DefinitionExemption implements Comparable<DefinitionExemption> {
  @Expose
  @SerializedName("gene")
  private final String m_gene;
  @Expose
  @SerializedName("requiredPositions")
  private final SortedSet<Long> m_requiredPositions;
  @Expose
  @SerializedName("ignoredPositions")
  private final SortedSet<VariantLocus> m_ignoredPositions;
  @Expose
  @SerializedName("extraPositions")
  private final SortedSet<VariantLocus> m_extraPositions;
  @Expose
  @SerializedName("ignoredAlleles")
  private final SortedSet<String> m_ignoredAlleles;
  @Expose
  @SerializedName("ignoredAllelesLc")
  private final SortedSet<String> m_ignoredAllelesLc;
  @Expose
  @SerializedName("unphasedDiplotypePriorities")
  private final SortedSet<UnphasedDiplotypePriority> m_unphasedDiplotypePriorities = new TreeSet<>();
  @Expose
  @SerializedName("amp1Alleles")
  private final List<String> m_amp1Alleles;
  @Expose
  @SerializedName("amp1Positions")
  private final SortedSet<Long> m_amp1Positions = new TreeSet<>();



  public DefinitionExemption(String gene, @Nullable SortedSet<Long> requiredPositions,
      @Nullable SortedSet<VariantLocus> ignoredPositions, @Nullable SortedSet<VariantLocus> extraPositions,
      @Nullable SortedSet<String> ignoredAlleles, @Nullable List<String> amp1Alleles) {
    m_gene = gene;
    m_requiredPositions = Objects.requireNonNullElse(requiredPositions, Collections.emptySortedSet());
    m_ignoredPositions = Objects.requireNonNullElse(ignoredPositions, Collections.emptySortedSet());
    m_extraPositions = Objects.requireNonNullElse(extraPositions, Collections.emptySortedSet());
    if (ignoredAlleles == null) {
      m_ignoredAlleles = Collections.emptySortedSet();
      m_ignoredAllelesLc = m_ignoredAlleles;
    } else {
      m_ignoredAlleles = ignoredAlleles;
      m_ignoredAllelesLc = ignoredAlleles.stream()
          .map(String::toLowerCase)
          .collect(Collectors.toCollection(TreeSet::new));
    }
    m_amp1Alleles = Objects.requireNonNullElse(amp1Alleles, Collections.emptyList());
  }


  public String getGene() {
    return m_gene;
  }


  /**
   * Gets the required positions that must be present before {@link NamedAlleleMatcher} will make a call.
   */
  public SortedSet<Long> getRequiredPositions() {
    return m_requiredPositions;
  }

  public boolean hasRequiredPositions() {
    return !m_requiredPositions.isEmpty();
  }

  public boolean isRequiredPosition(long pos) {
    return m_requiredPositions.contains(pos);
  }


  /**
   * Gets the positions from the original definition that should be ignored.
   * These get removed by the {@link org.pharmgkb.pharmcat.util.DataManager} when definitions are first pulled back from
   * PharmGKB.
   */
  public SortedSet<VariantLocus> getIgnoredPositions() {
    return m_ignoredPositions;
  }

  /**
   * Checks if the given position should be ignored.
   * <p>
   * <b>Currently only checks based on RSID!</b>
   */
  boolean shouldIgnorePosition(VariantLocus position) {
    return m_ignoredPositions.stream()
        .anyMatch(vl -> {
          if (vl.getRsid() != null) {
            return vl.getRsid().equals(position.getRsid());
          }
          return vl.getChromosomeHgvsName().equals(position.getChromosomeHgvsName());
        });
  }


  /**
   * Gets the extra positions to pull allele information for.
   */
  public SortedSet<VariantLocus> getExtraPositions() {
    return m_extraPositions;
  }


  /**
   * Gets the named alleles from the original definition that should be ignored.
   * These get removed by the {@link org.pharmgkb.pharmcat.util.DataManager} when definitions are first pulled back
   * from PharmGKB.
   */
  public SortedSet<String> getIgnoredAlleles() {
    return m_ignoredAlleles;
  }

  /**
   * Checks if the given named allele should be ignored.
   */
  boolean shouldIgnoreAllele(String allele) {
    return m_ignoredAllelesLc.contains(allele.toLowerCase());
  }


  public boolean hasUnphasedDiplotypePriorities() {
    return !m_unphasedDiplotypePriorities.isEmpty();
  }

  public SortedSet<UnphasedDiplotypePriority> getUnphasedDiplotypePriorities() {
    return m_unphasedDiplotypePriorities;
  }

  /**
   * Adds a priority diplotype when unphased data produces multiple calls.
   *
   * @param list a list of diplotypes; build using {@link #generateUnphasedPriorityKey(SortedSet)}
   * @param pick the diplotype to select given the {@code list} of diplotypes
   */
  public void addUnphasedDiplotypePriority(SortedSet<String> list, String pick) {
    UnphasedDiplotypePriority udp = new UnphasedDiplotypePriority(list, pick);
    for (UnphasedDiplotypePriority existing : m_unphasedDiplotypePriorities) {
      if (existing.getList().containsAll(udp.getList())) {
        throw new IllegalArgumentException("Cannot add " + String.join(", ", udp.getList()) +
            "- an existing superset of diplotypes already exist (" + String.join(", ", existing.getList()) + ")");
      }
      if (udp.getList().containsAll(existing.getList())) {
        throw new IllegalArgumentException("Cannot add " + String.join(", ", udp.getList()) +
            "- an existing subset of diplotypes already exist (" + String.join(", ", existing.getList()) + ")");
      }
    }
    m_unphasedDiplotypePriorities.add(udp);
  }


  public List<String> getAmp1Alleles() {
    return m_amp1Alleles;
  }

  public boolean hasAmp1Positions() {
    return m_amp1Positions.isEmpty();
  }

  public boolean isAmp1Position(long pos) {
    return m_amp1Positions.contains(pos);
  }

  public void addAmp1Position(VariantLocus pos) {
    m_amp1Positions.add(pos.getPosition());
  }


  @Override
  public int compareTo(DefinitionExemption o) {
    Preconditions.checkNotNull(m_gene);
    Preconditions.checkNotNull(o.getGene());
    return m_gene.compareTo(o.getGene());
  }


  public static String generateUnphasedPriorityKey(SortedSet<String> diplotypes) {
    return String.join("|", diplotypes);
  }

  public static String[] splitUnphasedPriorityKey(String key) {
    return key.split("\\|");
  }
}
