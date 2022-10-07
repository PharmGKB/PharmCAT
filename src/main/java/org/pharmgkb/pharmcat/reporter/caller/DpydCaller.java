package org.pharmgkb.pharmcat.reporter.caller;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.atomic.AtomicBoolean;
import com.google.common.collect.ImmutableList;
import org.apache.commons.lang3.ObjectUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Observation;


/**
 * This class handles calling DPYD diplotypes.
 *
 * @author Mark Woon
 */
public class DpydCaller {
  private static final String GENE = "DPYD";


  public static boolean isDpyd(String gene) {
    return GENE.equals(gene);
  }

  public static boolean isDpyd(GeneReport geneReport) {
    return GENE.equals(geneReport.getGene());
  }


  /**
   * Checks if DPYD was called with a true diplotype.
   */
  public static boolean hasTrueDiplotype(GeneCall geneCall) {
    if (!geneCall.isEffectivelyPhased()) {
      return false;
    }
    if (geneCall.getDiplotypes().size() > 1) {
      throw new IllegalStateException("Least function gene cannot have more than 1 diplotype");
    }
    return geneCall.getDiplotypes().size() == 1;
  }


  /**
   * Infer diplotypes from matcher results (based on true diplotypes).
   */
  public static List<Diplotype> inferFromDiplotypes(Collection<DiplotypeMatch> matches, Env env, DataSource source) {
    if (matches.size() == 0) {
      return ImmutableList.of(DiplotypeFactory.makeUnknownDiplotype(GENE, env, source));
    }
    List<Diplotype> diplotypes = new ArrayList<>();
    for (DiplotypeMatch dm : matches) {
      List<String> hapNames1 = new ArrayList<>(dm.getHaplotype1().getHaplotypeNames());
      List<String> hapNames2 = new ArrayList<>(dm.getHaplotype2().getHaplotypeNames());
      diplotypes.add(inferPhasedDiplotype(hapNames1, hapNames2, env, source));
    }
    return diplotypes;
  }

  /**
   * Infer diplotypes from matcher results (based on haplotype matches).
   */
  public static List<Diplotype> inferFromHaplotypeMatches(Collection<HaplotypeMatch> matches, Env env, DataSource source) {

    if (matches.size() == 0) {
      return ImmutableList.of(DiplotypeFactory.makeUnknownDiplotype(GENE, env, source));
    }
    List<String> hapNames = new ArrayList<>();
    for (HaplotypeMatch hm : matches) {
      hapNames.addAll(hm.getHaplotypeNames());
    }
    List<Diplotype> diplotypes = new ArrayList<>();
    diplotypes.add(inferUnphasedDiplotype(hapNames, env, source));
    return diplotypes;
  }

  /**
   * Infer diplotypes from outside call.
   */
  public static List<Diplotype> inferFromOutsideCall(String diplotype, Env env, DataSource source) {
    if (diplotype == null) {
      return ImmutableList.of(DiplotypeFactory.makeUnknownDiplotype(GENE, env, source));
    }

    String[] haplotypes = DiplotypeFactory.splitDiplotype(GENE, diplotype);
    List<String> hapNames1 = new ArrayList<>(DiplotypeFactory.splitHaplotype(haplotypes[0]));
    List<String> hapNames2 = new ArrayList<>();
    if (haplotypes.length == 2) {
      hapNames2.addAll(DiplotypeFactory.splitHaplotype(haplotypes[1]));
    }
    List<Diplotype> diplotypes = new ArrayList<>();
    diplotypes.add(inferPhasedDiplotype(hapNames1, hapNames2, env, source));
    return diplotypes;
  }

  private static Diplotype inferUnphasedDiplotype(List<String> hapNames, Env env, DataSource source) {
    Object[] hapData = makeHaplotypes(hapNames, env, source);
    //noinspection unchecked
    List<Haplotype> haplotypes = (List<Haplotype>)hapData[0];
    boolean isInferred = (Boolean)hapData[1];
    Haplotype hap1 = haplotypes.get(0);
    Haplotype hap2 = null;
    if (haplotypes.size() > 1) {
      hap2 = haplotypes.get(1);
    }
    Diplotype diplotype = new Diplotype(GENE, hap1, hap2);
    DiplotypeFactory.fillDiplotype(diplotype, env, source);
    if (isInferred || haplotypes.size() > 2) {
      diplotype.setObserved(Observation.INFERRED);
    }
    return diplotype;
  }

  private static Diplotype inferPhasedDiplotype(List<String> hapNames1, List<String> hapNames2, Env env,
      DataSource source) {

    Object[] hapData1 = makeHaplotypes(hapNames1, env, source);
    boolean isInferred = (Boolean)hapData1[1];
    //noinspection unchecked
    Haplotype hap1 = ((List<Haplotype>)hapData1[0]).get(0);
    Haplotype hap2 = null;
    if (hapNames2 != null && hapNames2.size() > 0) {
      Object[] hapData2 = makeHaplotypes(hapNames2, env, source);
      //noinspection unchecked
      hap2 = ((List<Haplotype>)hapData2[0]).get(0);
      isInferred = isInferred || (Boolean)hapData2[1];
    }
    Diplotype diplotype = new Diplotype(GENE, hap1, hap2);
    DiplotypeFactory.fillDiplotype(diplotype, env, source);
    if (isInferred || hapNames1.size() > 1 || (hapNames2 != null && hapNames2.size() > 1)) {
      diplotype.setObserved(Observation.INFERRED);
    }
    return diplotype;
  }


  private static Object[] makeHaplotypes(List<String> hapNames, Env env, DataSource source) {
    String refAllele = env.getReferenceAllele(GENE);
    AtomicBoolean inferred = new AtomicBoolean(false);
    List<Haplotype> haplotypes = hapNames.stream()
        .map(h -> {
          if (source == DataSource.DPWG) {
            GenePhenotype cpicGp = Objects.requireNonNull(env.getPhenotype(GENE, DataSource.CPIC));
            if ("normal function".equalsIgnoreCase(cpicGp.getHaplotypeFunction(h)) && !h.equals(refAllele)) {
              inferred.set(true);
              return env.makeHaplotype(GENE, refAllele, source);
            }
          }
          return env.makeHaplotype(GENE, h, source);
        })
        .sorted(new DpydActivityComparator(env))
        .toList();
    return new Object[] {
        haplotypes,
        inferred.get(),
    };
  }


  static class DpydActivityComparator implements Comparator<Haplotype> {
    private final Env m_env;

    public DpydActivityComparator(Env env) {
      m_env = env;
    }

    @Override
    public int compare(Haplotype o1, Haplotype o2) {
      if (o1 == o2) {
        return 0;
      }
      if (o1 == null) {
        return -1;
      } else if (o2 == null) {
        return 1;
      }
      int rez = ObjectUtils.compare(o1.getGene(), o2.getGene());
      if (rez != 0) {
        return rez;
      }

      // use activity scores from CPIC
      GenePhenotype cpicGp = Objects.requireNonNull(m_env.getPhenotype(GENE, DataSource.CPIC));
      rez = compare(cpicGp.getHaplotypeActivityScore(o1.getName()), cpicGp.getHaplotypeActivityScore(o2.getName()));
      if (rez != 0) {
        return rez;
      }

      // if same score, prefer one that's in DPWG
      GenePhenotype dpwgGp = Objects.requireNonNull(m_env.getPhenotype(GENE, DataSource.DPWG));
      int f1 = dpwgGp.getHaplotypes().containsKey(o1.getName()) ? 0 : 1;
      int f2 = dpwgGp.getHaplotypes().containsKey(o2.getName()) ? 0 : 1;
      rez = Integer.compare(f1, f2);
      if (rez != 0) {
        return rez;
      }

      return HaplotypeNameComparator.getComparator().compare(o1.getName(), o2.getName());
    }

    public static int compare(@Nullable Float a, @Nullable Float b) {
      if (a == null && b == null) {
        return 0;
      }
      if (a == null) {
        // b != null
        return -1;
      } else {
        // a != null
        if (b == null) {
          return 1;
        } else {
          // b != null
          return (a.compareTo(b));
        }
      }
    }
  }
}