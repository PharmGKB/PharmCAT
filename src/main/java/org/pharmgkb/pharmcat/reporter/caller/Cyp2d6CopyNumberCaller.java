package org.pharmgkb.pharmcat.reporter.caller;

import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.google.common.base.Preconditions;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Observation;


/**
 * This caller handles copy numbers for CYP2D6.
 *
 * @author Mark Woon
 */
public class Cyp2d6CopyNumberCaller {
  public static final String GENE = "CYP2D6";
  public static final String GTE = "≥";
  private static final Pattern sf_copyNumberPattern = Pattern.compile("\\*(\\d+)[Xx](≥?)(\\d+)");
  private static final Set<Integer> m_gteThree = new HashSet<>();


  public static void initialize(Env env) {
    if (m_gteThree.size() > 0) {
      return;
    }

    SortedSet<String> alleles = new TreeSet<>();
    alleles.addAll(Objects.requireNonNull(env.getPhenotype(Cyp2d6CopyNumberCaller.GENE, DataSource.CPIC))
        .getHaplotypes().keySet());
    alleles.addAll(Objects.requireNonNull(env.getPhenotype(Cyp2d6CopyNumberCaller.GENE, DataSource.DPWG))
        .getHaplotypes().keySet());

    for (String allele : alleles) {
      if (allele.contains(GTE)) {
        Matcher m = sf_copyNumberPattern.matcher(allele);
        if (!m.matches()) {
          throw new UnsupportedOperationException(allele + " does not match expected pattern");
        }
        if (!m.group(3).equals("3")) {
          throw new UnsupportedOperationException("Expecting x" + GTE + "3 but got " + m.group(2) + m.group(3));
        }
        m_gteThree.add(Integer.parseInt(m.group(1)));
      }
    }
  }


  public static Diplotype inferDiplotype(GeneReport report, Diplotype diplotype, Env env, DataSource source) {
    Preconditions.checkArgument(diplotype.getGene().equals(GENE), "Can only be used on CYP2D6");

    if (!report.isOutsideCall()) {
      return diplotype;
    }

    if (diplotype.isPhenotypeOnly() || diplotype.isOutsideActivityScore()) {
      return diplotype;
    }

    Object[] r1 = inferHaplotype(diplotype.getAllele1());
    boolean needsInfer1 = (Boolean)r1[0];
    Object[] r2 = inferHaplotype(diplotype.getAllele2());
    boolean needsInfer2 = (Boolean)r2[0];

    if (!needsInfer1 && !needsInfer2) {
      return diplotype;
    }

    Haplotype hap1 = needsInfer1 ? env.makeHaplotype(GENE, (String)r1[1], source) : diplotype.getAllele1();
    Haplotype hap2 = needsInfer2 ? env.makeHaplotype(GENE, (String)r2[1], source) : diplotype.getAllele2();
    Diplotype inferredDiplotype = new Diplotype(GENE, hap1, hap2, env, source);
    inferredDiplotype.setObserved(Observation.INFERRED);
    return inferredDiplotype;
  }


  private static Object[] inferHaplotype(Haplotype haplotype) {
    if (haplotype == null) {
      return new Object[] {false, null};
    }
    Matcher m = sf_copyNumberPattern.matcher(haplotype.getName());
    if (!m.matches()) {
      return new Object[] {false, haplotype.getName()};
    }
    Integer hap = Integer.parseInt(m.group(1));
    if (!m_gteThree.contains(hap)) {
      return new Object[] {false, haplotype.getName()};
    }
    int cn = Integer.parseInt(m.group(3));
    if (cn <= 2) {
      return new Object[] {false, haplotype.getName()};
    }
    if (haplotype.getName().contains(GTE)) {
      if (cn == 3) {
        // ignore >= 3
        return new Object[] {false, haplotype.getName()};
      }
    }
    return new Object[] {true, "*" + hap + "x" + GTE + "3"};
  }
}
