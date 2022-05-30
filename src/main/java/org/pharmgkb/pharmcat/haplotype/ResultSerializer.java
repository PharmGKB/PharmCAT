package org.pharmgkb.pharmcat.haplotype;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.apache.commons.io.IOUtils;
import org.apache.commons.text.StringSubstitutor;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.haplotype.model.Variant;


/**
 * Serializer/Deserializer for {@link Result}.
 *
 * @author Mark Woon
 */
public class ResultSerializer {
  private static final Gson sf_gson = new GsonBuilder().serializeNulls().excludeFieldsWithoutExposeAnnotation()
      .setPrettyPrinting().create();
  private boolean m_alwaysShowUnmatchedHaplotypes;
  private final SimpleDateFormat m_dateFormat = new SimpleDateFormat("MM/dd/yy");
  private String m_htmlTemplate;


  public ResultSerializer() {
  }


  public ResultSerializer alwaysShowUnmatchedHaplotypes(boolean alwaysShowUnmatchedHaplotypes) {
    m_alwaysShowUnmatchedHaplotypes = alwaysShowUnmatchedHaplotypes;
    return this;
  }


  private String getHtmlTemplate() throws IOException {
    if (m_htmlTemplate == null) {
      try (BufferedReader reader = Files.newBufferedReader(PathUtils.getPathToResource(getClass(), "template.html"))) {
        m_htmlTemplate = IOUtils.toString(reader);
      }
    }
    return m_htmlTemplate;
  }


  public ResultSerializer toJson(Result result, Path jsonFile) throws IOException {
    Preconditions.checkNotNull(result);
    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"), "Output JSON file needs to end in '.json'");

    try (BufferedWriter writer = Files.newBufferedWriter(jsonFile, StandardCharsets.UTF_8)) {
      writer.write(sf_gson.toJson(result));
    }
    return this;
  }


  public Result fromJson(Path jsonFile) throws IOException {
    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"));
    Preconditions.checkArgument(Files.isRegularFile(jsonFile));

    try (BufferedReader reader = Files.newBufferedReader(jsonFile, StandardCharsets.UTF_8)) {
      return sf_gson.fromJson(reader, Result.class);
    }
  }



  public ResultSerializer toHtml(Result result, Path htmlFile) throws IOException {
    Preconditions.checkNotNull(result);
    Preconditions.checkNotNull(htmlFile);
    Preconditions.checkArgument(htmlFile.toString().endsWith(".html"));

    StringBuilder builder = new StringBuilder();
    for (GeneCall call : result.getGeneCalls()) {
      MatchData matchData = call.getMatchData();
      Map<Long, String> refAlleleMap = new HashMap<>();
      Arrays.stream(matchData.getPositions()).forEach(vl -> refAlleleMap.put(vl.getPosition(), vl.getRef()));
      Set<Long> highlightPositions = new HashSet<>();

      builder.append("<h3>")
          .append(call.getGene())
          .append("</h3>");

      builder.append("<ul>");
      for (DiplotypeMatch diplotype : call.getDiplotypes()) {
        builder.append("<li>")
            .append(diplotype.getName())
            .append(" (")
            .append(diplotype.getScore())
            .append(")</li>");
      }
      builder.append("</ul>");

      builder.append("<table class=\"table table-striped table-hover table-condensed\">");
      // position
      builder.append("<tr>");
      builder.append("<th>Definition Position</th>");
      for (Variant v : call.getVariants()) {
        builder.append("<th>")
            .append(v.getPosition())
            .append("</th>");
      }
      builder.append("</tr>");
      // rsid
      builder.append("<tr>");
      builder.append("<th></th>");
      for (Variant v : call.getVariants()) {
        builder.append("<th>");
        if (v.getRsid() != null) {
          builder.append(v.getRsid());
        }
        builder.append("</th>");
      }
      builder.append("</tr>");
      // VCF position
      builder.append("<tr>");
      builder.append("<th>VCF Position</th>");
      for (Variant v : call.getVariants()) {
        builder.append("<th>")
            .append(v.getPosition())
            .append("</th>");
      }
      builder.append("</tr>");
      // sample
      builder.append("<tr>");
      builder.append("<th>VCF REF,ALTs</th>");
      for (Variant v : call.getVariants()) {
        builder.append("<th>")
            .append(v.getVcfAlleles())
            .append("</th>");
      }
      builder.append("</tr>");

      builder.append("<tr class=\"success\">");
      builder.append("<th>VCF Call</th>");
      for (Variant v : call.getVariants()) {
        if (v.getVcfCall() != null) {
          SortedSet<String> alleles = new TreeSet<>(Splitter.on(v.isPhased() ? "|" : "/").splitToList(v.getVcfCall()));
          boolean isNonRef = false;
          if (alleles.size() > 1 || !alleles.first().equals(refAlleleMap.get(v.getPosition()))) {
            isNonRef = true;
            highlightPositions.add(v.getPosition());
          }
          builder.append("<th");
          if (isNonRef) {
            builder.append(" class=\"text-danger\"");
          }
          builder.append(">")
              .append(v.getVcfCall())
              .append("</th>");
        } else {
          builder.append("<th />");
        }
      }
      builder.append("</tr>");

      Set<String> matchedHaplotypeNames = new HashSet<>();
      if (call.getHaplotypes().size() > 0) {
        for (BaseMatch hm : call.getHaplotypes()) {
          matchedHaplotypeNames.add(hm.getHaplotype().getName());
          printAllele(builder, hm.getHaplotype().getName(), hm.getHaplotype().getPermutations().pattern(), "info",
              highlightPositions);
          for (String seq : hm.getSequences()) {
            printAllele(builder, null, seq, null, highlightPositions);
          }
        }
      }
      if (m_alwaysShowUnmatchedHaplotypes || matchedHaplotypeNames.size() == 0) {
        for (NamedAllele haplotype : matchData.getHaplotypes()) {
          if (!matchedHaplotypeNames.contains(haplotype.getName())) {
            printAllele(builder, haplotype.getName(), haplotype.getPermutations().pattern(), "danger", highlightPositions);
          }
        }
      }

      builder.append("</table>");

      if (matchData.getMissingPositions().size() > 0) {
        builder.append("<p>There ");
        if (matchData.getMissingPositions().size() > 1) {
          builder.append("were ");
        } else {
          builder.append("was ");
        }
        builder.append(matchData.getMissingPositions().size())
            .append(" missing positions from the VCF file:</p>")
            .append("<ul>");
        for (VariantLocus variant : matchData.getMissingPositions()) {
          builder.append("<li>")
              .append(variant.getPosition())
              .append(" (")
              .append(variant.getChromosomeHgvsName())
              .append(")</li>");
        }
        builder.append("</ul>");

        if (call.getUncallableHaplotypes().size() > 0) {
          builder.append("<p>The following haplotype(s) were eliminated from consideration:</p>")
              .append("<ul>");
          for (String name : call.getUncallableHaplotypes()) {
            builder.append("<li>")
                .append(name)
                .append("</li>");
          }
          builder.append("</ul>");
        }

        if (call.getHaplotypes().size() > 0) {
          builder.append("<p>The following haplotypes were called even though tag positions were missing:</p>")
              .append("<ul>");
          for (BaseMatch hm : call.getHaplotypes()) {
            if (hm.getHaplotype().getMissingPositions().size() > 0) {
              builder.append("<li>Called ")
                  .append(hm.getName())
                  .append(" without ")
                  .append(hm.getHaplotype().getMissingPositions().stream()
                      .map(VariantLocus::getChromosomeHgvsName)
                      .collect(Collectors.joining(", ")))
                  .append("</li>");
            }
          }
          builder.append("</ul>");
        }
      }
    }

    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(htmlFile, StandardCharsets.UTF_8))) {
      Map<String, String> varMap = new HashMap<>();
      varMap.put("title", "PharmCAT Allele Call Report for " + result.getMetadata().getInputFilename());
      varMap.put("content", builder.toString());
      varMap.put("timestamp", m_dateFormat.format(new Date()));
      StringSubstitutor sub = new StringSubstitutor(varMap);
      writer.println(sub.replace(getHtmlTemplate()));
    }
    return this;
  }


  private void printAllele(StringBuilder builder, @Nullable String name, String allele,
      @Nullable String rowClass, Set<Long> highlightPositions) {

    SortedSet<Variant> variants = new TreeSet<>();
    for (String posAllele : allele.split(";")) {
      String[] parts = posAllele.split(":");
      String a = parts[1];
      if (a.equals(".?")) {
        a = "";
      }
      long vcfPosition = Long.parseLong(parts[0]);
      variants.add(new Variant(vcfPosition, null, a, ""));
    }

    builder.append("<tr");
    if (rowClass != null) {
      builder.append(" class=\"")
          .append(rowClass)
          .append("\"");
    }
    builder.append("><th>");
    if (name != null) {
      builder.append(name);
    }
    builder.append("</th>");

    for (Variant variant : variants) {
      String vcfCall = variant.getVcfCall();
      if (vcfCall != null && vcfCall.contains("\\")) {
        vcfCall = vcfCall.replaceAll("\\\\", "");
      }
      builder.append("<td");
      boolean isAny = ".*?".equals(vcfCall);
      if (highlightPositions.contains(variant.getPosition()) && !isAny) {
        builder.append(" class=\"text-danger\"");
      }
      builder.append(">");
      if (name == null || isAny) {
        builder.append(vcfCall);
      } else {
        builder.append("<b>")
            .append(vcfCall)
            .append("</b>");
      }
      builder.append("</td>");
    }

    builder.append("</tr>");
  }
}
