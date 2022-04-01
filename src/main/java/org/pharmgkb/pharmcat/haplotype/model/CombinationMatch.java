package org.pharmgkb.pharmcat.haplotype.model;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Splitter;
import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;


/**
 * This class represents a combination {@link NamedAllele} match.
 *
 * @author Mark Woon
 */
public class CombinationMatch extends BaseMatch {
  public static final String COMBINATION_JOINER = " + ";
  public static final Splitter COMBINATION_NAME_SPLITTER = Splitter.on(COMBINATION_JOINER).trimResults();
  private final SortedSet<NamedAllele> m_componentHaplotypes = new TreeSet<>();
  private final VariantLocus[] m_refVariants;


  public CombinationMatch(VariantLocus[] refVariants, NamedAllele namedAllele, String sequence) {
    m_refVariants = refVariants;
    m_componentHaplotypes.add(namedAllele);
    setName(buildName());
    setHaplotype(buildHaplotype());
    addSequence(sequence);
  }

  public CombinationMatch(CombinationMatch combinationMatch) {
    m_refVariants = combinationMatch.getRefVariants();
    m_componentHaplotypes.addAll(combinationMatch.getComponentHaplotypes());
    setName(buildName());
    setHaplotype(buildHaplotype());
    addSequence(combinationMatch.getSequences().first());
  }


  private String buildName() {
    return m_componentHaplotypes.stream()
        .map(NamedAllele::getName)
        .collect(Collectors.joining(COMBINATION_JOINER));
  }

  private NamedAllele buildHaplotype() {
    if (m_componentHaplotypes.first().getAlleles().length != m_componentHaplotypes.first().getCpicAlleles().length) {
      throw new IllegalStateException(m_componentHaplotypes.first() + " has different number of alleles and cpicAlleles");
    }
    String id = m_componentHaplotypes.stream().map(NamedAllele::getId).collect(Collectors.joining(COMBINATION_JOINER));
    String[] alleles = new String[m_componentHaplotypes.first().getAlleles().length];
    String[] cpicAlleles = new String[alleles.length];
    for (int x = 0; x < alleles.length; x += 1) {
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
    NamedAllele na = new NamedAllele(id, getName(), alleles, cpicAlleles, false);
    na.initialize(m_refVariants);
    na.setCombinationOrPartial(m_componentHaplotypes.size() > 1);
    return na;
  }


  SortedSet<NamedAllele> getComponentHaplotypes() {
    return m_componentHaplotypes;
  }

  private VariantLocus[] getRefVariants() {
    return m_refVariants;
  }


  public boolean canMerge(NamedAllele namedAllele, String seq) {
    if (getSequences().first().equals(seq)) {
      for (NamedAllele na : m_componentHaplotypes) {
        for (int x = 0; x < namedAllele.getAlleles().length; x += 1) {
          // haplotypes cannot have overlapping variations (no sharing!)
          if (namedAllele.getAlleles()[x] != null && na.getAlleles()[x] != null) {
            return false;
          }
        }
      }
      return true;
    }
    return false;
  }

  public void merge(NamedAllele namedAllele) {
    m_componentHaplotypes.add(namedAllele);
    setName(buildName());
    setHaplotype(buildHaplotype());
  }


  @Override
  public int compareTo(BaseMatch o) {
    if (this == o) {
      return 0;
    }
    if (o instanceof CombinationMatch cm) {
      Iterator<NamedAllele> it1 = m_componentHaplotypes.iterator();
      Iterator<NamedAllele> it2 = cm.getComponentHaplotypes().iterator();
      while (it1.hasNext()) {
        NamedAllele na1 = it1.next();
        if (it2.hasNext()) {
          NamedAllele na2 = it2.next();
          int rez = ObjectUtils.compare(na1, na2);
          if (rez != 0) {
            return rez;
          }
        } else {
          return 1;
        }
      }
      if (it2.hasNext()) {
        return -1;
      }
      return compareSequences(o);
    }
    if (o instanceof HaplotypeMatch hm) {
      int rez = ObjectUtils.compare(m_componentHaplotypes.first(), hm.getHaplotype());
      if (rez != 0) {
        return rez;
      }
      return 1;
    }
    return ObjectUtils.compare(getName(), o.getName());
  }

  @Override
  public int hashCode() {
    return Objects.hash(getName(), Arrays.hashCode(getSequences().toArray()));
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (!(obj instanceof CombinationMatch that)) {
      return false;
    }
    return Objects.equals(getName(), that.getName()) &&
        Arrays.equals(getSequences().toArray(), that.getSequences().toArray());
  }
}
