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
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;


/**
 * This caller handles copy numbers for CYP2D6.
 *
 * @author Mark Woon
 */
public class Cyp2d6CopyNumberCaller {
  public static final String GENE = "CYP2D6";
  private static final Pattern sf_copyNumberPattern = Pattern.compile("\\*(\\d+)[Xx](" + TextConstants.GTE + "?)(\\d+)");
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
      if (allele.contains(TextConstants.GTE)) {
        Matcher m = sf_copyNumberPattern.matcher(allele);
        if (!m.matches()) {
          throw new UnsupportedOperationException(allele + " does not match expected pattern");
        }
        if (!m.group(3).equals("3")) {
          throw new UnsupportedOperationException("Expecting x" + TextConstants.GTE + "3 but got " + m.group(2) +
              m.group(3));
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
    inferredDiplotype.setInferred(true);
    inferredDiplotype.setInferredSourceDiplotype(diplotype);
    return inferredDiplotype;
  }


  private static Object[] inferHaplotype(Haplotype haplotype) {
    if (haplotype == null) {
      return new Object[] {false, null};
    }
    String name = inferHaplotypeName(haplotype.getName());
    return new Object[] {!name.equals(haplotype.getName()), name};
  }


  public static String inferHaplotypeName(String haplotypeName) {

    Matcher m = sf_copyNumberPattern.matcher(haplotypeName);
    if (!m.matches()) {
      return haplotypeName;
    }
    Integer hap = Integer.parseInt(m.group(1));
    if (!m_gteThree.contains(hap)) {
      return haplotypeName;
    }
    int cn = Integer.parseInt(m.group(3));
    if (cn <= 2) {
      return haplotypeName;
    }
    if (haplotypeName.contains(TextConstants.GTE)) {
      if (cn == 3) {
        // ignore >= 3
        return haplotypeName;
      }
    }
    return "*" + hap + "x" + TextConstants.GTE + "3";
  }
}
