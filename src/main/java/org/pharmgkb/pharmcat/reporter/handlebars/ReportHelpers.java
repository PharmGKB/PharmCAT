package org.pharmgkb.pharmcat.reporter.handlebars;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
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
      if (builder.length() > 0) {
        builder.append(", ");
      }
      builder.append("<a href=\"")
          .append(url)
          .append("\">source");
      if (urls.size() > 1) {
          builder.append(" ")
              .append(x);
      }
      builder.append("</a>");
    }
    return builder.toString();
  }


  public static String printRecMap(Map<String, String> data) {
    if (data.size() == 1) {
      String val = data.values().iterator().next();
      if (val.equalsIgnoreCase(TextConstants.NA)) {
        val = TextConstants.NA;
      }
      return "<p>" + val + "</p>";
    }

    StringBuilder builder = new StringBuilder()
        .append("<dl class=\"compact mt-0\">");
    for (String key : data.keySet()) {
      builder.append("<dt>")
          .append(key)
          .append(":</dt><dd>");
      String val = data.get(key);
      if (val.equalsIgnoreCase(TextConstants.NA)) {
        val = TextConstants.NA;
      }
      builder.append(val);
      builder.append("</dd>");
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
      return diplotype.printFunctionPhrase();
    }
  }

  public static String gsPhenotype(Diplotype diplotype) {
    return diplotype.printPhenotypes();
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


  public static String amdSubtitle(GeneReport geneReport) {
    StringBuilder builder = new StringBuilder();

    if (isDpyd(geneReport.getGene()) && geneReport.getComponentDiplotypes().size() == 0) {
      builder.append("Haplotype");
    } else {
      builder.append("Genotype");
    }
    if (geneReport.getMatcherDiplotypes().size() > 1) {
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

  public static List<String> amdMessages(GeneReport geneReport) {
    return geneReport.getMessages().stream()
        .filter(MessageAnnotation.isMessage)
        .map(MessageAnnotation::getMessage)
        .collect(Collectors.toList());
  }

  public static List<String> amdExtraPositionNotes(GeneReport geneReport) {
    return geneReport.getMessages().stream()
        .filter(MessageAnnotation.isExtraPositionNote)
        .map(MessageAnnotation::getMessage)
        .collect(Collectors.toList());
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
}
