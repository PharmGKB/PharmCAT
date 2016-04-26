package org.pharmgkb.pharmcat.haplotype;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.text.StrSubstitutor;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.json.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.json.HaplotyperResult;
import org.pharmgkb.pharmcat.haplotype.model.json.Metadata;


/**
 * @author Mark Woon
 */
public class Report {
  private DefinitionReader m_definitionReader;
  private Path m_vcfFile;
  private HaplotyperResult m_root = new HaplotyperResult();


  public Report(@Nonnull DefinitionReader definitionReader) {
    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
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


  protected Report gene(@Nonnull String gene, @Nonnull List<DiplotypeMatch> matches,
      Collection<SampleAllele> sampleAlleles) {

    Preconditions.checkNotNull(gene);
    Preconditions.checkNotNull(matches);

    GeneCall geneCall = new GeneCall();
    geneCall.setGene(gene);

    // get haplotype/diplotype info
    for (DiplotypeMatch dm : matches) {
      geneCall.addDiplotype(dm);
    }

    // get position info
    Map<Integer, SampleAllele> alleleMap = new HashMap<>();
    for (SampleAllele allele : sampleAlleles) {
      alleleMap.put(allele.getPosition(), allele);
    }
    for (Variant variant : m_definitionReader.getHaplotypePositions().get(gene)) {
      SampleAllele allele = alleleMap.get(variant.getPosition());
      String call;
      if (allele.isPhased()) {
        call = allele.getAllele1() + "|" + allele.getAllele2();
      } else {
        call = allele.getAllele1() + "/" + allele.getAllele2();
      }
      geneCall.add(new org.pharmgkb.pharmcat.haplotype.model.json.Variant(variant.getPosition(), variant.getRsid(), call));
    }

    //geneCall.setHaplotypesNotCalled();
    DefinitionFile tsvFile = m_definitionReader.getDefinitionFiles().get(gene);
    geneCall.setGeneVersion(tsvFile.getContentVersion() + " (" + tsvFile.getContentDate() + ")");
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
      for (org.pharmgkb.pharmcat.haplotype.model.json.Variant v : call.getVariants()) {
        builder.append("<th>")
            .append(v.getPosition())
            .append("</th>");
      }
      builder.append("</tr>");
      // rsid
      builder.append("<tr>");
      builder.append("<th></th>");
      for (org.pharmgkb.pharmcat.haplotype.model.json.Variant v : call.getVariants()) {
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
      for (org.pharmgkb.pharmcat.haplotype.model.json.Variant v : call.getVariants()) {
        builder.append("<th>")
            .append(v.getVcfCall())
            .append("</th>");
      }
      builder.append("</tr>");

      for (HaplotypeMatch hm : call.getHaplotypes()) {
        printAllele(builder, hm.getHaplotype().getName(), hm.getHaplotype().getPermutations().pattern());
        for (String seq : hm.getSequences()) {
          printAllele(builder, null, seq);
        }
        builder.append("</tr>");
      }
      builder.append("</table>");

    }

    System.out.println("Printing to " + htmlFile);
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(htmlFile, StandardCharsets.UTF_8))) {
      Map<String, String> varMap = new HashMap<>();
      varMap.put("title", "PharmCAT Allele Call Report for " + m_root.getMetadata().getInputFile());
      varMap.put("content", builder.toString());
      varMap.put("timestamp", SimpleDateFormat.getDateTimeInstance().format(new Date()));
      StrSubstitutor sub = new StrSubstitutor(varMap);
      String template = IOUtils.toString(getClass().getResourceAsStream("template.html"));
      writer.println(sub.replace(template));
    }
  }


  private void printAllele(StringBuilder builder, String name, String allele) {

    if (name != null) {
      builder.append("<tr class=\"info\">")
          .append("<th>")
          .append(name);
    } else {
      builder.append("<tr>")
          .append("<th>");
    }
    builder.append("</th>");

    for (String part : allele.split(";")) {
      String a = part.split(":")[1];
      if (a.equals(".?")) {
        a = "";
      }
      builder.append("<td>");
      if (name == null) {
        builder.append(a);
      } else {
        builder.append("<b>")
            .append(a)
            .append("</b>");
      }
      builder.append("</td>");
    }

    builder.append("</tr>");
  }
}
