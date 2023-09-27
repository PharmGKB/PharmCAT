package org.pharmgkb.pharmcat.reporter.caller;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicBoolean;
import org.apache.commons.lang3.ObjectUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * This class handles calling DPYD diplotypes.
 *
 * @author Mark Woon
 */
public class DpydCaller {
  public static final String GENE = "DPYD";


  public static boolean isDpyd(String gene) {
    return GENE.equals(gene);
  }

  public static boolean isDpyd(GeneReport geneReport) {
    return GENE.equals(geneReport.getGene());
  }


  /**
   * Infer diplotypes from matcher results (based on true diplotypes).
   */
  public static List<Diplotype> inferFromDiplotypes(Collection<DiplotypeMatch> matches, Env env, DataSource source) {
    List<Diplotype> diplotypes = new ArrayList<>();
    if (matches.isEmpty()) {
      diplotypes.add(DiplotypeFactory.makeUnknownDiplotype(GENE, env, source));
    } else {
      for (DiplotypeMatch dm : matches) {
        List<String> hapNames1 = new ArrayList<>(dm.getHaplotype1().getHaplotypeNames());
        List<String> hapNames2;
        if (dm.getHaplotype2() == null) {
          hapNames2 = Collections.emptyList();
        } else {
          hapNames2 = new ArrayList<>(dm.getHaplotype2().getHaplotypeNames());
        }
        diplotypes.add(inferPhasedDiplotype(hapNames1, hapNames2, env, source));
      }
    }

    DiplotypeFactory diplotypeFactory = new DiplotypeFactory(GENE, env);
    SortedSet<Diplotype> srcDiplotypes = new TreeSet<>(diplotypeFactory.makeDiplotypes(matches, source));
    diplotypes.forEach(d -> d.setInferredSourceDiplotypes(srcDiplotypes));
    return diplotypes;
  }

  /**
   * Infer diplotypes from matcher results (based on haplotype matches).
   */
  public static List<Diplotype> inferFromHaplotypeMatches(Collection<HaplotypeMatch> matches, Env env, DataSource source) {

    List<Diplotype> diplotypes = new ArrayList<>();
    if (matches.isEmpty()) {
      diplotypes.add(DiplotypeFactory.makeUnknownDiplotype(GENE, env, source));
    } else {
      List<String> hapNames = new ArrayList<>();
      for (HaplotypeMatch hm : matches) {
        hapNames.addAll(hm.getHaplotypeNames());
      }
      diplotypes.add(inferUnphasedDiplotype(hapNames, env, source));
    }

    DiplotypeFactory diplotypeFactory = new DiplotypeFactory(GENE, env);
    SortedSet<Diplotype> srcDiplotypes = new TreeSet<>(diplotypeFactory.makeDiplotypesFromHaplotypeMatches(matches, source));
    diplotypes.forEach(d -> d.setInferredSourceDiplotypes(srcDiplotypes));
    return diplotypes;
  }

  /**
   * Infer diplotypes from outside call.
   */
  public static List<Diplotype> inferFromOutsideCall(OutsideCall outsideCall, Env env, DataSource source) {
    String diplotype = outsideCall.getDiplotype();
    Diplotype inferredDiplotype;
    if (diplotype == null) {
      inferredDiplotype = DiplotypeFactory.makeUnknownDiplotype(GENE, env, source);
    } else {
      String[] haplotypes = DiplotypeFactory.splitDiplotype(GENE, diplotype);
      Arrays.sort(haplotypes, HaplotypeNameComparator.getComparator());
      List<String> hapNames1 = new ArrayList<>(DiplotypeFactory.splitHaplotype(haplotypes[0]));
      List<String> hapNames2 = new ArrayList<>();
      if (haplotypes.length == 2) {
        hapNames2.addAll(DiplotypeFactory.splitHaplotype(haplotypes[1]));
      }
      inferredDiplotype = inferPhasedDiplotype(hapNames1, hapNames2, env, source);
    }

    inferredDiplotype.setInferredSourceDiplotype(new Diplotype(outsideCall, env, source));
    return List.of(inferredDiplotype);
  }

  private static Diplotype inferUnphasedDiplotype(List<String> hapNames, Env env, DataSource source) {
    InferredHaps hapData = makeHaplotypes(hapNames, env, source);
    List<Haplotype> haplotypes = hapData.haplotypes;
    boolean isInferred = hapData.isInferred;
    Haplotype hap1 = haplotypes.get(0);
    Haplotype hap2 = null;
    if (haplotypes.size() > 1) {
      hap2 = haplotypes.get(1);
    }
    Diplotype diplotype = new Diplotype(GENE, hap1, hap2, env, source);
    if (isInferred || haplotypes.size() > 2) {
      diplotype.setInferred(true);
    }
    return diplotype;
  }

  private static Diplotype inferPhasedDiplotype(List<String> hapNames1, List<String> hapNames2, Env env,
      DataSource source) {

    InferredHaps hapData1 = makeHaplotypes(hapNames1, env, source);
    boolean isInferred = hapData1.isInferred;
    Haplotype hap1 = hapData1.haplotypes.get(0);
    Haplotype hap2 = null;
    if (hapNames2 != null && !hapNames2.isEmpty()) {
      InferredHaps hapData2 = makeHaplotypes(hapNames2, env, source);
      hap2 = hapData2.haplotypes.get(0);
      isInferred = isInferred || hapData2.isInferred;
    }
    Diplotype diplotype = new Diplotype(GENE, hap1, hap2, env, source);
    if (isInferred || hapNames1.size() > 1 || (hapNames2 != null && hapNames2.size() > 1)) {
      diplotype.setInferred(true);
    }
    return diplotype;
  }


  private static InferredHaps makeHaplotypes(List<String> hapNames, Env env, DataSource source) {
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
    return new InferredHaps(haplotypes, inferred.get());
  }


  /**
   * Wrapper class to hold multiple returned values.
   */
  private static class InferredHaps {
    List<Haplotype> haplotypes;
    boolean isInferred;

    private InferredHaps(List<Haplotype> haplotypes, boolean isInferred) {
      this.haplotypes = haplotypes;
      this.isInferred = isInferred;
    }
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
