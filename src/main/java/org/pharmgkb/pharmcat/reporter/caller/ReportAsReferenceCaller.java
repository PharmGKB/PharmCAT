package org.pharmgkb.pharmcat.reporter.caller;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.util.DataManager;
import org.pharmgkb.pharmcat.util.DataSerializer;


/**
 * This caller handles the report-as-reference mechanic.
 *
 * @author Mark Woon
 */
public class ReportAsReferenceCaller {
  public static final String DATA_FILE_NAME = "reportAsReference.json";

  private static Map<String, List<String>> m_reportAsReference;


  private static @Nullable List<String> getReportAsReference(String gene) {
    if (m_reportAsReference == null) {
      Path file = DataManager.DEFAULT_DEFINITION_DIR.resolve(DATA_FILE_NAME);
      System.out.println("Reading report-as-reference data in " + file);

      try (BufferedReader reader = Files.newBufferedReader(file)) {
        //noinspection unchecked
        m_reportAsReference = (Map<String, List<String>>)DataSerializer.GSON.fromJson(reader, Map.class);
      } catch (IOException ex) {
        throw new IllegalStateException("Error reading " + file, ex);
      }
    }
    return m_reportAsReference.get(gene);
  }

  private static boolean shouldOverride(String gene) {
    return getReportAsReference(gene) != null;
  }

  private static boolean shouldOverride(String gene, String namedAllele) {
    List<String> rez = getReportAsReference(gene);
    if (rez == null) {
      return false;
    }
    return rez.contains(namedAllele);
  }


  public static Set<DiplotypeMatch> override(Env env, String gene, Set<DiplotypeMatch> matches) {

    if (!shouldOverride(gene)) {
      return matches;
    }

    LinkedHashSet<DiplotypeMatch> updated = new LinkedHashSet<>();
    for (DiplotypeMatch dm : matches) {
      if (shouldOverride(gene, dm.getHaplotype1()) ||
          shouldOverride(gene, dm.getHaplotype2())) {
        updated.add(new DiplotypeMatch(
            shouldOverride(gene, dm.getHaplotype1())
                ? override(env, gene, dm.getHaplotype1())
                : dm.getHaplotype1(),
            shouldOverride(gene, dm.getHaplotype2())
                ? override(env, gene, dm.getHaplotype2())
                : dm.getHaplotype2(),
            dm.getDataset()));
      } else {
        updated.add(dm);
      }
    }
    return updated;
  }

  public static List<HaplotypeMatch> override(Env env, String gene, List<HaplotypeMatch> matches) {

    if (!shouldOverride(gene)) {
      return matches;
    }

    List<HaplotypeMatch> updated = new ArrayList<>();
    for (HaplotypeMatch hm : matches) {
      if (shouldOverride(gene, hm)) {
        updated.add(override(env, gene, hm));
      } else {
        updated.add(hm);
      }
    }
    return updated;
  }


  public static OutsideCall override(Env env, OutsideCall oc) {

    if (!shouldOverride(oc.getGene())) {
      return oc;
    }

    String diplotype = oc.getHaplotypes().stream()
        .map(h -> {
          if (shouldOverride(oc.getGene(), h)) {
            return env.getDefinitionReader().getDefinitionFile(oc.getGene()).getReferenceNamedAllele().getName();
          }
          return h;
        })
        .collect(Collectors.joining("/ "));
    StringBuilder builder = new StringBuilder()
        .append(oc.getGene())
        .append("\t")
        .append(diplotype);
    if (oc.getPhenotype() != null) {
      builder.append(oc.getPhenotype());
    }
    builder.append("\t");
    if (oc.getActivityScore() != null) {
      builder.append(oc.getActivityScore());
    }
    return new OutsideCall(env, builder.toString(), 999);
  }


  private static boolean shouldOverride(String gene, BaseMatch bm) {

    if (bm instanceof HaplotypeMatch hm) {
      return shouldOverride(gene, hm.getName());

    } else if (bm instanceof CombinationMatch cm) {
      for (NamedAllele na : cm.getComponentHaplotypes()) {
        if (shouldOverride(gene, na.getName())) {
          return true;
        }
      }
      return false;

    }
    throw new IllegalStateException("What is this: " + bm.getClass().getName());
  }


  private static <T extends BaseMatch> T override(Env env, String gene, T bm) {

    if (bm instanceof HaplotypeMatch hm) {
      HaplotypeMatch match = new HaplotypeMatch(env.getDefinitionReader().getReferenceHaplotype(gene));
      hm.getSequences().forEach(match::addSequence);
      //noinspection unchecked
      return (T)match;

    } else if (bm instanceof CombinationMatch cm) {
      CombinationMatch match = new CombinationMatch(cm.getRefVariants(),
          env.getDefinitionReader().getReferenceHaplotype(gene), cm.getSequences().first());
      for (NamedAllele na : cm.getComponentHaplotypes()) {
        if (!shouldOverride(gene, na.getName())) {
          match.merge(na);
        }
      }
      //noinspection unchecked
      return (T)match;
    }
    throw new IllegalStateException("What is this: " + bm.getClass().getName());
  }


  public static Collection<Diplotype> override(Env env, DataSource source, String gene, Collection<Diplotype> matches) {

    if (!shouldOverride(gene)) {
      return matches;
    }

    LinkedHashSet<Diplotype> diplotypes = new LinkedHashSet<>();
    for (Diplotype dip : matches) {
      if (shouldOverride(dip.getAllele1()) ||
          shouldOverride(dip.getAllele2())) {
        diplotypes.add(new Diplotype(gene,
            override(dip.getAllele1()),
            override(dip.getAllele2()),
            env, source));
      } else {
        diplotypes.add(dip);
      }
    }
    return diplotypes;
  }

  public static Diplotype override(Env env, DataSource source, Diplotype dip) {

    if (!shouldOverride(dip.getGene())) {
      return dip;
    }
    if (shouldOverride(dip.getAllele1()) ||
        shouldOverride(dip.getAllele2())) {
      return new Diplotype(dip.getGene(),
          override(dip.getAllele1()),
          override(dip.getAllele2()),
          env, source);
    }
    return dip;
  }


  private static boolean shouldOverride(Haplotype h1) {
    if (h1 == null) {
      return false;
    }
    return shouldOverride(h1.getGene(), h1.getName());
  }

  private static @Nullable String override(Haplotype h1) {
    if (h1 == null) {
      return null;
    }
    if (shouldOverride(h1.getGene(), h1.getName())) {
      return "*1";
    }
    return h1.getName();
  }


  public static Optional<MessageAnnotation> addMessage(GeneReport geneReport, Collection<?> orig,
      Collection<?> current) {
    return addMessage(geneReport,
        orig.stream().map(o -> o.toString().replaceAll(geneReport.getGene() + ":", ""))
            .collect(Collectors.joining(", ")),
        current.stream().map(o -> o.toString().replaceAll(geneReport.getGene() + ":", ""))
            .collect(Collectors.joining(", ")));
  }

  public static Optional<MessageAnnotation> addMessage(GeneReport geneReport, String orig, String current) {
    if (orig.equals(current)) {
      return Optional.empty();
    }
    MessageAnnotation msg = new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "report-as-reference",
        "Rewriting " + orig + " as " + current);
    geneReport.addMessage(msg);
    return Optional.of(msg);
  }
}
