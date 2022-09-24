package org.pharmgkb.pharmcat.reporter.handlebars;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.cpic.Publication;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationReport;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;

import static org.pharmgkb.pharmcat.reporter.DpydCaller.isDpyd;
import static org.pharmgkb.pharmcat.reporter.TextConstants.UNCALLED;


/**
 * Class to hold methods that can be used as helper methods in the Handlebars templating system.
 *
 * @author Ryan Whaley
 */
@SuppressWarnings("unused")
public class ReportHelpers {
  private static final String sf_variantAlleleTemplate = "<td class=\"%s\">%s%s</td>";

  public static String listGenes(Collection<String> genes) {
    StringBuilder builder = new StringBuilder();
    for (String gene : genes) {
      if (builder.length() > 0) {
        builder.append(", ");
      }
      builder.append("<a href=\"#")
          .append(gene)
          .append("\">")
          .append(gene)
          .append("</a>");
    }
    return builder.toString();
  }

  public static String listSources(Collection<String> urls) {
    StringBuilder builder = new StringBuilder();
    int x = 0;
    for (String url : urls) {
      x += 1;
      if (x > 1) {
        builder.append(", ");
      }
      builder.append(externalHref(url, Integer.toString(x)));
    }
    return builder.toString();
  }


  public static String printRecMap(Map<String, String> data) {
    if (data.size() == 1) {
      return "<p>" + capitalizeNA(data.values().iterator().next()) + "</p>";
    }

    StringBuilder builder = new StringBuilder()
        .append("<dl class=\"compact mt-0\">");
    for (String key : data.keySet()) {
      builder.append("<dt>")
          .append(key)
          .append(":</dt><dd>")
          .append(capitalizeNA(data.get(key)))
          .append("</dd>");
    }
    builder.append("</dl>");
    return builder.toString();
  }


  public static String pluralize(String word, Object col) {
    if (moreThanOne(col)) {
      return word + "s";
    }
    return word;
  }

  public static String variantAlleles(VariantReport variantReport) {
    String cellStyle = variantReport.isNonwildtype() ? "nonwild" : "";
    String mismatch = variantReport.isMismatch() ? "<div class=\"callMessage\">Mismatch: Called allele does not match allele definitions</div>" : "";
    if (variantReport.isMismatch()) {
      cellStyle = StringUtils.strip(cellStyle + " mismatch");
    }
    return String.format(sf_variantAlleleTemplate, cellStyle, variantReport.getCall(), mismatch);
  }


  public static boolean gsShowDrug(String drug, Collection<String> drugs) {
    return drugs.contains(drug);
  }


  public static String gsCall(Diplotype diplotype) {
    if (diplotype.isUnknownAlleles()) {
      return UNCALLED;
    }
    return diplotype.printDisplay();
  }

  public static String gsFunction(Diplotype diplotype) {
    if (diplotype.isCombination() && isDpyd(diplotype.getGene())) {
      return TextConstants.SEE_DRUG;
    } else {
      return capitalizeNA(diplotype.printFunctionPhrase());
    }
  }

  public static String gsPhenotype(Diplotype diplotype) {
    return capitalizeNA(diplotype.printPhenotypes());
  }


  public static String rxAnnotationClass(DataSource source, String drug) {
    StringBuilder builder = new StringBuilder();
    if (source == DataSource.CPIC) {
      builder.append("cpic-");
    } else {
      builder.append("dpwg-");
    }
    builder.append(drug);
    return builder.toString();
  }

  public static boolean rxIsCpicWarfarin(String drug, DataSource source) {
    return drug.equals("warfarin") && source == DataSource.CPIC;
  }

  public static boolean rxDpydInferred(Genotype genotype) {
    return genotype.isInferred() && genotype.getDiplotypes().stream()
        .map(Diplotype::getGene)
        .anyMatch(g -> g.equals("DPYD"));
  }

  public static boolean rxInferred(Genotype genotype) {
    return genotype.isInferred() && genotype.getDiplotypes().stream()
        .map(Diplotype::getGene)
        .noneMatch(g -> g.equals("DPYD"));
  }

  public static String rxGenotype(Genotype genotype, AnnotationReport annotationReport) {
    if (genotype.getDiplotypes().size() == 0) {
      return TextConstants.UNKNOWN_GENOTYPE;
    }
    StringBuilder builder = new StringBuilder();
    for (Diplotype diplotype : genotype.getDiplotypes()) {
      if (builder.length() > 0) {
        builder.append(";<br />");
      }
      builder
          .append("<a href=\"#")
          .append(diplotype.getGene())
          .append("\">")
          .append(diplotype.getGene())
          .append("</a>:");
      String call = diplotype.printBare();
      if (call.length() <= 15) {
        builder.append(call);
      } else {
        int idx = call.indexOf("/");
        if (idx == -1) {
          builder.append("<br />")
              .append(call);
        } else {
          String a = call.substring(0, idx + 1);
          String b = call.substring(idx + 1);
          if (a.length() > 15) {
            builder.append("<br />");
          }
          builder.append(a)
              .append("<br />")
              .append(b);
        }
      }
    }
    if (annotationReport.getHighlightedVariants().size() > 0) {
      for (String var : annotationReport.getHighlightedVariants()) {
        if (builder.length() > 0) {
          builder.append(";<br />");
        }
        builder.append("<span id=\"")
            .append(annotationReport.getLocalId())
            .append("-")
            .append(var, 0, var.indexOf(":"))
            .append("\">")
            .append(var)
            .append("</span>");
      }
    }
    return builder.toString();
  }

  public static List<MessageAnnotation> rxAnnotationMessages(AnnotationReport annotationReport) {
    return annotationReport.getMessages().stream()
        .filter(MessageAnnotation.isMessage)
        .toList();
  }


  public static String amdSubtitle(GeneReport geneReport) {
    StringBuilder builder = new StringBuilder();

    if (isDpyd(geneReport.getGene()) && geneReport.getMatcherComponentDiplotypes().size() == 0) {
      builder.append("Haplotype");
    } else {
      builder.append("Genotype");
    }
    if (geneReport.getSourceDiplotypes().size() > 1) {
      builder.append("s");
    }
    builder.append(" ");
    if (geneReport.isOutsideCall()) {
      builder.append("Reported");
    } else {
      builder.append("Matched");
    }
    return builder.toString();
  }

  public static boolean amdNoCall(GeneReport report) {
    if (report.getCallSource() == CallSource.NONE) {
      return true;
    } else if (report.getCallSource() == CallSource.MATCHER) {
      if (report.getVariantReports().size() == 0 ||
          report.getVariantReports().stream().allMatch(VariantReport::isMissing)) {
        return true;
      }
    }
    return false;
  }

  public static boolean amdIsSingleCall(GeneReport report) {
    return report.printDisplayCalls().size() == 1;
  }

  public static List<String> amdGeneCalls(GeneReport report) {
    return report.printDisplayCalls();
  }

  public static String amdGeneCall(GeneReport report) {
    return report.printDisplayCalls().get(0);
  }


  public static String amdPhaseStatus(GeneReport geneReport) {
    if (geneReport.isOutsideCall()) {
      return "Unavailable for calls made outside PharmCAT";
    }
    return geneReport.isPhased() ? "Phased" : "Unphased";
  }

  public static boolean amdShowUnphasedNote(GeneReport geneReport) {
    return !geneReport.isPhased() && !isDpyd(geneReport.getGeneDisplay());
  }

  public static boolean amdHasUncalledHaps(GeneReport geneReport) {
    return geneReport.getUncalledHaplotypes() != null &&
        geneReport.getUncalledHaplotypes().size() > 0;
  }

  public static String amdUncalledHaps(GeneReport geneReport) {
    return String.join(", ", geneReport.getUncalledHaplotypes());
  }


  public static long amdTotalMissingVariants(GeneReport geneReport) {
    return geneReport.getVariantReports().stream()
        .filter(VariantReport::isMissing)
        .count();
  }

  public static int amdTotalVariants(GeneReport geneReport) {
    return geneReport.getVariantReports().size();
  }

  public static List<MessageAnnotation> amdMessages(GeneReport geneReport) {
    return geneReport.getMessages().stream()
        .filter(MessageAnnotation.isMessage)
        .toList();
  }

  public static List<String> amdExtraPositionNotes(GeneReport geneReport) {
    return geneReport.getMessages().stream()
        .filter(MessageAnnotation.isExtraPositionNote)
        .map(MessageAnnotation::getMessage)
        .toList();
  }


  public static boolean moreThanOne(Object obj) {
    if (obj == null) {
      return false;
    }
    //noinspection rawtypes
    if (obj instanceof Collection col) {
      return col.size() > 1;
    }
    //noinspection rawtypes
    if (obj instanceof Map map) {
      return map.size() > 1;
    }
    return false;
  }

  public static boolean thereAre(Object obj) {
    if (obj == null) {
      return false;
    }
    //noinspection rawtypes
    if (obj instanceof Collection col) {
      return col.size() > 0;
    }
    //noinspection rawtypes
    if (obj instanceof Map map) {
      return map.size() > 0;
    }
    return false;
  }

  public static int add(int x, int y) {
    return x + y;
  }

  public static String externalHref(String href, String text) {
    return "<a href=\"" + href + "\" target=\"_blank\" rel=\"noopener noreferrer\">" + text + "</a>";
  }

  private static final Pattern sf_endsWithPuncuationPattern = Pattern.compile(".*?\\p{Punct}$");

  public static String printCitation(Publication pub) {
    String url = "https://www.ncbi.nlm.nih.gov/pubmed/" + pub.getPmid();
    String period = "";
    if (!sf_endsWithPuncuationPattern.matcher(pub.getTitle()).matches()) {
      period = ".";
    }
    return externalHref(url, pub.getTitle()) + period + " <i>" + pub.getJournal() + "</i>. " + pub.getYear() +
        ". PMID:" + pub.getPmid();
  }

  public static String capitalizeNA(String text) {
    if (TextConstants.NA.equalsIgnoreCase(text)) {
      return "N/A";
    }
    return text;
  }

  public static String messageClass(MessageAnnotation msg) {
    return msg.getName()
        .replaceAll("\\p{Punct}", "-")
        .replaceAll("\\s+", "-")
        .replaceAll("-+", "-");
  }
}
