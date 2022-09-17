package org.pharmgkb.pharmcat.reporter;

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
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
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


  public static List<Diplotype> inferDiplotypes(@SuppressWarnings("rawtypes") Collection matches,
      boolean hasTrueDiplotype, Env env, DataSource source) {

    if (matches.size() > 0) {
      List<Diplotype> diplotypes = new ArrayList<>();
      if (hasTrueDiplotype) {
        for (Object obj : matches) {
          List<String> hapNames1 = new ArrayList<>();
          List<String> hapNames2 = new ArrayList<>();
          if (obj instanceof DiplotypeMatch dm) {
            // from matcher
            hapNames1.addAll(dm.getHaplotype1().getHaplotypeNames());
            hapNames2.addAll(dm.getHaplotype2().getHaplotypeNames());
          } else if (obj instanceof String text) {
            // from outside call
            String[] haplotypes = DiplotypeFactory.splitDiplotype(GENE, text);
            hapNames1.addAll(DiplotypeFactory.splitHaplotype(haplotypes[0]));
            if (haplotypes.length == 2) {
              hapNames2.addAll(DiplotypeFactory.splitHaplotype(haplotypes[1]));
            }
          } else {
            throw new IllegalStateException("Unexpected type: " + obj.getClass());
          }
          diplotypes.add(inferPhasedDiplotype(hapNames1, hapNames2, env, source));
        }
      } else {
        List<String> hapNames = new ArrayList<>();
        for (Object obj : matches) {
          if (obj instanceof BaseMatch bm) {
            // from matcher
            hapNames.addAll(bm.getHaplotypeNames());
          } else {
            throw new IllegalStateException("Unexpected type: " + obj.getClass());
          }
        }
        diplotypes.add(inferUnphasedDiplotype(hapNames, env, source));
      }
      return diplotypes;
    } else {
      return ImmutableList.of(DiplotypeFactory.makeUnknownDiplotype(GENE, env, source));
    }
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


  static Object[] makeHaplotypes(List<String> hapNames, Env env, DataSource source) {
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
