package org.pharmgkb.pharmcat.haplotype.model;

import java.util.Collection;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.CombinationMatcher;


/**
 * This class represents a combination {@link NamedAllele} match.
 *
 * @author Mark Woon
 */
public class CombinationMatch extends BaseMatch {
  @Expose
  @SerializedName("componentHaplotypes")
  private final SortedSet<NamedAllele> m_componentHaplotypes = new TreeSet<>();
  @Expose
  @SerializedName("partials")
  private final SortedMap<Long, String> m_partials = new TreeMap<>();
  private final VariantLocus[] m_refVariants;


  /**
   * Constructor for creating a copy of an existing {@link CombinationMatch}.
   */
  public CombinationMatch(CombinationMatch combinationMatch) {
    m_refVariants = combinationMatch.getRefVariants();
    m_componentHaplotypes.addAll(combinationMatch.getComponentHaplotypes());
    setName(buildName());
    setHaplotype(buildHaplotype(false));
    addSequence(combinationMatch.getSequences().first());
  }

  /**
   * Primary constructor.
   */
  public CombinationMatch(VariantLocus[] refVariants, String seq, Collection<NamedAllele> components,
      @Nullable Map<Long, String> partials) {
    m_refVariants = refVariants;
    addSequence(seq);
    m_componentHaplotypes.addAll(components);
    if (partials != null) {
      partials.keySet().forEach(p -> m_partials.put(p, partials.get(p)));
    }
    setName(buildName());
    setHaplotype(buildHaplotype(false));
  }


  public int getNumCombinations() {
    return m_componentHaplotypes.size();
  }

  private String buildName() {
    StringBuilder builder = new StringBuilder();
    int count = 0;
    for (NamedAllele na : m_componentHaplotypes) {
      if (na.isReference()) {
        continue;
      }
      if (!builder.isEmpty()) {
        builder.append(CombinationMatcher.COMBINATION_JOINER);
      }
      if (na.isCombination()) {
        builder.append(CombinationMatcher.extractCombinationName(na.getName()));
        count += 2;
      } else {
        builder.append(na.getName());
        count += 1;
      }
    }
    for (long pos : m_partials.keySet()) {
      if (!builder.isEmpty()) {
        builder.append(CombinationMatcher.COMBINATION_JOINER);
      }
      builder.append(m_partials.get(pos));
      count += 1;
    }
    if (count > 1) {
      return "[" + builder + "]";
    }
    return builder.toString();
  }

  /**
   * Builds a new {@link NamedAllele} based on component haplotypes.
   *
   * @param isOffReferencePartial if true, will set the score to 0
   */
  private NamedAllele buildHaplotype(boolean isOffReferencePartial) {
    int length = m_refVariants.length;
    for (NamedAllele na : m_componentHaplotypes) {
      if (na.getAlleles().length != length) {
        throw new IllegalStateException("Component haplotypes have unexpected number of alleles");
      }
      if (na.getCpicAlleles().length != length) {
        throw new IllegalStateException(na + " has different number of alleles and cpicAlleles");
      }
    }
    StringBuilder idBuilder = new StringBuilder();
    SortedSet<VariantLocus> missingPositions = new TreeSet<>();
    for (NamedAllele na : m_componentHaplotypes) {
      if (!idBuilder.isEmpty()) {
        idBuilder.append(CombinationMatcher.COMBINATION_JOINER);
      }
      idBuilder.append(na.getId());
      missingPositions.addAll(na.getMissingPositions());
    }
    String[] alleles = new String[length];
    String[] cpicAlleles = new String[length];
    for (int x = 0; x < length; x += 1) {
      for (NamedAllele na : m_componentHaplotypes) {
        if (alleles[x] == null) {
          alleles[x] = na.getAlleles()[x];
        } else if (na.getAlleles()[x] != null && !alleles[x].equals(na.getAlleles()[x])) {
          throw new IllegalStateException(getName() + " has different alleles @ index " + x);
        }
        if (cpicAlleles[x] == null) {
          cpicAlleles[x] = na.getCpicAlleles()[x];
        } else if (na.getCpicAlleles()[x] != null && !cpicAlleles[x].equals(na.getCpicAlleles()[x])) {
          throw new IllegalStateException(getName() + " has different CPIC alleles @ index " + x);
        }
      }
    }
    NamedAllele na = new NamedAllele(idBuilder.toString(), getName(), alleles, cpicAlleles, missingPositions,
        false, m_componentHaplotypes.size(), m_partials.size());
    if (isOffReferencePartial) {
      na.initialize(m_refVariants, 0);
    } else {
      na.initialize(m_refVariants);
    }
    return na;
  }


  public SortedSet<NamedAllele> getComponentHaplotypes() {
    return m_componentHaplotypes;
  }

  private VariantLocus[] getRefVariants() {
    return m_refVariants;
  }
}
