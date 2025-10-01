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
import java.util.stream.Stream;
import org.apache.commons.lang3.ObjectUtils;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * This class handles calling lowest-function gene diplotypes.
 *
 * @author Mark Woon
 */
public class LowestFunctionGeneCaller {
  public static final String DPYD = "DPYD";


  /**
   * Infer diplotypes from matcher results (based on true diplotypes).
   */
  public static List<Diplotype> inferFromDiplotypes(String gene, Env env, DataSource source,
      DiplotypeFactory diplotypeFactory, Collection<DiplotypeMatch> matches) {

    if (matches.isEmpty()) {
      return List.of(DiplotypeFactory.makeUnknownDiplotype(gene, env, source));
    }

    SortedSet<Diplotype> diplotypes;
    if (gene.equals("DPYD")) {
      diplotypes = new TreeSet<>(new DiplotypeCpicActivityScoreComparator(env, gene));
    } else if (gene.equals("RYR1")) {
      diplotypes = new TreeSet<>(Ryr1DiplotypeMalignancyComparator.INSTANCE);
    } else {
      throw new IllegalArgumentException("Don't know how to deal with " + gene);
    }
    for (DiplotypeMatch dm : matches) {
      List<String> hapNames1 = new ArrayList<>(dm.getHaplotype1().getHaplotypeNames());
      List<String> hapNames2;
      if (dm.getHaplotype2() == null) {
        hapNames2 = Collections.emptyList();
      } else {
        hapNames2 = new ArrayList<>(dm.getHaplotype2().getHaplotypeNames());
      }
      diplotypes.add(inferPhasedDiplotype(gene, hapNames1, hapNames2, env, source));
    }

    Diplotype finalDiplotype = diplotypes.first();
    finalDiplotype.setInferredSourceDiplotypes(new TreeSet<>(diplotypeFactory.makeDiplotypes(matches, source)));
    return List.of(finalDiplotype);
  }

  /**
   * Infer diplotypes from matcher results (based on haplotype matches).
   */
  public static List<Diplotype> inferFromHaplotypeMatches(String gene, Env env, DataSource source,
      DiplotypeFactory diplotypeFactory, Collection<HaplotypeMatch> matches) {

    List<Diplotype> diplotypes = new ArrayList<>();
    if (matches.isEmpty()) {
      diplotypes.add(DiplotypeFactory.makeUnknownDiplotype(gene, env, source));
    } else {
      List<String> hapNames = new ArrayList<>();
      for (HaplotypeMatch hm : matches) {
        hapNames.addAll(hm.getHaplotypeNames());
      }
      diplotypes.add(inferUnphasedDiplotype(gene, hapNames, env, source));
    }

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
      inferredDiplotype = DiplotypeFactory.makeUnknownDiplotype(outsideCall.getGene(), env, source);
    } else {
      String[] haplotypes = DiplotypeFactory.splitDiplotype(outsideCall.getGene(), diplotype);
      Arrays.sort(haplotypes, HaplotypeNameComparator.getComparator());
      List<String> hapNames1 = new ArrayList<>(DiplotypeFactory.splitHaplotype(haplotypes[0]));
      List<String> hapNames2 = new ArrayList<>();
      if (haplotypes.length == 2) {
        hapNames2.addAll(DiplotypeFactory.splitHaplotype(haplotypes[1]));
      }
      inferredDiplotype = inferPhasedDiplotype(outsideCall.getGene(), hapNames1, hapNames2, env, source);
    }

    inferredDiplotype.setInferredSourceDiplotype(new Diplotype(outsideCall, env, source));
    return List.of(inferredDiplotype);
  }

  private static Diplotype inferUnphasedDiplotype(String gene, List<String> hapNames, Env env, DataSource source) {
    InferredHaps hapData = makeHaplotypes(gene, hapNames, env, source);
    List<Haplotype> haplotypes = hapData.haplotypes;
    boolean isInferred = hapData.isInferred;
    Haplotype hap1 = haplotypes.get(0);
    Haplotype hap2 = null;
    if (haplotypes.size() > 1) {
      hap2 = haplotypes.get(1);
    }
    Diplotype diplotype = new Diplotype(gene, hap1, hap2, env, source);
    if (isInferred || haplotypes.size() > 2) {
      diplotype.setInferred(true);
    }
    return diplotype;
  }

  private static Diplotype inferPhasedDiplotype(String gene, List<String> hapNames1, List<String> hapNames2, Env env,
      DataSource source) {

    InferredHaps hapData1 = makeHaplotypes(gene, hapNames1, env, source);
    boolean isInferred = hapData1.isInferred;
    Haplotype hap1 = hapData1.haplotypes.get(0);
    Haplotype hap2 = null;
    if (!hapNames2.isEmpty()) {
      InferredHaps hapData2 = makeHaplotypes(gene, hapNames2, env, source);
      hap2 = hapData2.haplotypes.get(0);
      isInferred = isInferred || hapData2.isInferred;
    }
    Diplotype diplotype = new Diplotype(gene, hap1, hap2, env, source);
    if (isInferred || hapNames1.size() > 1 || hapNames2.size() > 1) {
      diplotype.setInferred(true);
    }
    return diplotype;
  }


  private static InferredHaps makeHaplotypes(String gene, List<String> hapNames, Env env, DataSource source) {
    String refAllele = env.getReferenceAllele(gene);
    AtomicBoolean inferred = new AtomicBoolean(false);
    Stream<Haplotype> hapStream = hapNames.stream()
        .map(h -> {
          if (source == DataSource.DPWG) {
            GenePhenotype cpicGp = Objects.requireNonNull(env.getPhenotype(gene, DataSource.CPIC));
            if ("normal function".equalsIgnoreCase(cpicGp.getHaplotypeFunction(h)) && !h.equals(refAllele)) {
              inferred.set(true);
              return env.makeHaplotype(gene, refAllele, source);
            }
          }
          return env.makeHaplotype(gene, h, source);
        });
    if (gene.equals("DPYD")) {
      hapStream = hapStream.sorted(new DpydActivityComparator(env));
    } else if (gene.equals("RYR1")) {
      hapStream = hapStream.sorted(Ryr1ActivityComparator.INSTANCE);
    }
    List<Haplotype> haplotypes = hapStream.toList();
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


  /**
   * Diplotype comparator that compares based on CPIC activity score.
   * Expects to sort the score ascending, so this gives haplotypes with no activity score a higher score so that it
   * comes later.
   */
  private static class DiplotypeCpicActivityScoreComparator implements Comparator<Diplotype> {
    private final GenePhenotype m_cpicGp;

    public DiplotypeCpicActivityScoreComparator(Env env, String gene) {
      m_cpicGp = Objects.requireNonNull(env.getPhenotype(gene, DataSource.CPIC));
    }

    private float getHaplotypeActivityScore(String hapName) {
      Float score = m_cpicGp.getHaplotypeActivityScore(hapName);
      // use large number to rank it low
      return Objects.requireNonNullElse(score, 2.0f);
    }

    @Override
    public int compare(@Nullable Diplotype o1, @Nullable Diplotype o2) {
      if (o1 == o2) {
        return 0;
      }
      if (o1 == null) {
        return -1;
      } else if (o2 == null) {
        return 1;
      }

      // use activity scores from CPIC
      float as1 = getHaplotypeActivityScore(Objects.requireNonNull(o1.getAllele1()).getName()) +
          getHaplotypeActivityScore(Objects.requireNonNull(o1.getAllele2()).getName());
      float as2 = getHaplotypeActivityScore(Objects.requireNonNull(o2.getAllele1()).getName()) +
          getHaplotypeActivityScore(Objects.requireNonNull(o2.getAllele2()).getName());

      int rez = ObjectUtils.compare(as1, as2);
      if (rez != 0) {
        return rez;
      }

      rez = HaplotypeNameComparator.getComparator().compare(o1.getAllele1().getName(), o2.getAllele1().getName());
      if (rez != 0) {
        return rez;
      }
      return HaplotypeNameComparator.getComparator().compare(o1.getAllele2().getName(), o2.getAllele2().getName());
    }
  }

  static class DpydActivityComparator implements Comparator<Haplotype> {
    private final GenePhenotype m_cpicGp;
    private final GenePhenotype m_dpwgGp;

    public DpydActivityComparator(Env env) {
      m_cpicGp = Objects.requireNonNull(env.getPhenotype("DPYD", DataSource.CPIC));
      m_dpwgGp = Objects.requireNonNull(env.getPhenotype("DPYD", DataSource.DPWG));
    }

    @Override
    public int compare(@Nullable Haplotype o1, @Nullable Haplotype o2) {
      if (o1 == o2) {
        return 0;
      }
      if (o1 == null) {
        return -1;
      } else if (o2 == null) {
        return 1;
      }

      // use activity scores from CPIC
      int rez = ObjectUtils.compare(m_cpicGp.getHaplotypeActivityScore(o1.getName()),
          m_cpicGp.getHaplotypeActivityScore(o2.getName()));
      if (rez != 0) {
        return rez;
      }

      // if same score, prefer one that's in DPWG
      int f1 = m_dpwgGp.getHaplotypes().containsKey(o1.getName()) ? 0 : 1;
      int f2 = m_dpwgGp.getHaplotypes().containsKey(o2.getName()) ? 0 : 1;
      rez = Integer.compare(f1, f2);
      if (rez != 0) {
        return rez;
      }

      return HaplotypeNameComparator.getComparator().compare(o1.getName(), o2.getName());
    }
  }


  private static class Ryr1DiplotypeMalignancyComparator implements Comparator<Diplotype> {
    static Ryr1DiplotypeMalignancyComparator INSTANCE = new Ryr1DiplotypeMalignancyComparator();

    @Override
    public int compare(@Nullable Diplotype o1, @Nullable Diplotype o2) {
      if (o1 == o2) {
        return 0;
      }
      if (o1 == null) {
        return -1;
      } else if (o2 == null) {
        return 1;
      }

      int m1 = checkMalignant(Objects.requireNonNull(o1.getAllele1())) + checkMalignant(Objects.requireNonNull(o1.getAllele2()));
      int m2 = checkMalignant(Objects.requireNonNull(o2.getAllele1())) + checkMalignant(Objects.requireNonNull(o2.getAllele2()));
      // rank more malignant first
      int rez = Integer.compare(m2, m1);
      if (rez != 0) {
        return rez;
      }

      rez = HaplotypeNameComparator.getComparator().compare(o1.getAllele1().getName(), o2.getAllele1().getName());
      if (rez != 0) {
        return rez;
      }
      return HaplotypeNameComparator.getComparator().compare(o1.getAllele2().getName(), o2.getAllele2().getName());
    }

    private int checkMalignant(Haplotype haplotype) {
      return haplotype.getFunction().contains("Malignant") ? 1 : 0;
    }
  }

  static class Ryr1ActivityComparator implements Comparator<Haplotype> {
    static Ryr1ActivityComparator INSTANCE = new Ryr1ActivityComparator();

    @Override
    public int compare(@Nullable Haplotype o1, @Nullable Haplotype o2) {
      if (o1 == o2) {
        return 0;
      }
      if (o1 == null) {
        return -1;
      } else if (o2 == null) {
        return 1;
      }

      boolean o1Malignant = o1.getFunction().contains("Malignant");
      boolean o2Malignant = o2.getFunction().contains("Malignant");
      if (o1Malignant != o2Malignant) {
        if (o1Malignant) {
          return -1;
        }
        return 1;
      }
      return HaplotypeNameComparator.getComparator().compare(o1.getName(), o2.getName());
    }
  }
}
