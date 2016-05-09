package org.pharmgkb.pharmcat.haplotype;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.text.StrSubstitutor;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.json.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.json.HaplotyperResult;
import org.pharmgkb.pharmcat.haplotype.model.json.Metadata;
import org.pharmgkb.pharmcat.haplotype.model.json.Variant;


/**
 * @author Mark Woon
 */
public class Report {
  private DefinitionReader m_definitionReader;
  private Path m_vcfFile;
  private HaplotyperResult m_root = new HaplotyperResult();
  private SimpleDateFormat m_dateFormat = new SimpleDateFormat("MM/dd/yy");
  private Map<String, MatchData> m_dataMap = new HashMap<>();


  public Report(@Nonnull DefinitionReader definitionReader) {
    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
  }

  HaplotyperResult getResults() {
    return m_root;
  }


  public Report forFile(@Nonnull Path vcfFile) {
    Preconditions.checkNotNull(vcfFile);
    Preconditions.checkArgument(vcfFile.toString().endsWith(".vcf"));
    Preconditions.checkArgument(Files.isRegularFile(vcfFile));

    m_vcfFile = vcfFile;
    Metadata md = new Metadata();
    //md.setCpicAnnotatorBuild();
    //md.setCpicDataBuild();
    md.setDate(new Date());
    //md.setGenomeAssembly();
    md.setInputFile(vcfFile.getName(vcfFile.getNameCount() - 1).toString());
    m_root.setMetadata(md);

    return this;
  }


  protected Report gene(@Nonnull String gene, @Nonnull List<DiplotypeMatch> matches, @Nonnull MatchData dataset) {

    Preconditions.checkNotNull(gene);
    Preconditions.checkNotNull(matches);

    m_dataMap.put(gene, dataset);

    GeneCall geneCall = new GeneCall();
    geneCall.setGene(gene);

    // get haplotype/diplotype info
    for (DiplotypeMatch dm : matches) {
      geneCall.addDiplotype(dm);
    }

    // get position info
    Map<Integer, SampleAllele> alleleMap = new HashMap<>();
    for (SampleAllele allele : dataset.geneSampleMap.values()) {
      alleleMap.put(allele.getPosition(), allele);
    }
    for (VariantLocus variant : dataset.positions) {
      SampleAllele allele = alleleMap.get(variant.getPosition());
      String call;
      if (allele.isPhased()) {
        call = allele.getAllele1() + "|" + allele.getAllele2();
      } else {
        call = allele.getAllele1() + "/" + allele.getAllele2();
      }
      geneCall.add(new Variant(variant.getPosition(), variant.getRsid(), call));
    }

    //geneCall.setHaplotypesNotCalled();
    DefinitionFile tsvFile = m_definitionReader.getDefinitionFile(gene);
    geneCall.setGeneVersion(tsvFile.getContentVersion() + " (" + m_dateFormat.format(tsvFile.getModificationDate()) + ")");
    geneCall.setChromosome(tsvFile.getChromosome());
    m_root.addDiplotypeCall(geneCall);

    return this;
  }


  public void print() throws IOException {

    Preconditions.checkState(m_vcfFile != null);
    Path jsonFile = m_vcfFile.getParent().resolve(PathUtils.getBaseFilename(m_vcfFile) + ".json");
    Gson gson = new GsonBuilder().setPrettyPrinting().create();

    try (BufferedWriter writer = Files.newBufferedWriter(jsonFile, StandardCharsets.UTF_8)) {
      writer.write(gson.toJson(m_root));
    }
  }


  public void printHtml() throws IOException {

    Preconditions.checkState(m_vcfFile != null);
    Path htmlFile = m_vcfFile.getParent().resolve(PathUtils.getBaseFilename(m_vcfFile) + ".html");

    StringBuilder builder = new StringBuilder();
    for (GeneCall call : m_root.getGeneCalls()) {
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
      builder.append("<th></th>");
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
      // sample
      builder.append("<tr class=\"success\">");
      builder.append("<th>VCF</th>");
      for (Variant v : call.getVariants()) {
        builder.append("<th>")
            .append(v.getVcfCall())
            .append("</th>");
      }
      builder.append("</tr>");

      if (call.getHaplotypes().size() > 0) {
        for (HaplotypeMatch hm : call.getHaplotypes()) {
          printAllele(builder, hm.getHaplotype().getName(), hm.getHaplotype().getPermutations().pattern(), "info");
          for (String seq : hm.getSequences()) {
            printAllele(builder, null, seq, null);
          }
        }
      } else {
        for (NamedAllele haplotype : m_dataMap.get(call.getGene()).haplotypes) {
          printAllele(builder, haplotype.getName(), haplotype.getPermutations().pattern(), "danger");
        }
      }

      builder.append("</table>");

    }

    System.out.println("Printing to " + htmlFile);
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(htmlFile, StandardCharsets.UTF_8))) {
      Map<String, String> varMap = new HashMap<>();
      varMap.put("title", "PharmCAT Allele Call Report for " + m_root.getMetadata().getInputFile());
      varMap.put("content", builder.toString());
      varMap.put("timestamp", m_dateFormat.format(new Date()));
      StrSubstitutor sub = new StrSubstitutor(varMap);
      String template = IOUtils.toString(getClass().getResourceAsStream("template.html"));
      writer.println(sub.replace(template));
    }
  }


  private void printAllele(@Nonnull StringBuilder builder, @Nullable String name, @Nonnull String allele,
      @Nullable String rowClass) {

    SortedSet<Variant> variants = new TreeSet<>();
    for (String posAllele : allele.split(";")) {
      String[] parts = posAllele.split(":");
      String a = parts[1];
      if (a.equals(".?")) {
        a = "";
      }
      variants.add(new Variant(Integer.parseInt(parts[0]), null, a));
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
      builder.append("<td>");
      if (name == null) {
        builder.append(variant.getVcfCall());
      } else {
        builder.append("<b>")
            .append(variant.getVcfCall())
            .append("</b>");
      }
      builder.append("</td>");
    }

    builder.append("</tr>");
  }
}
