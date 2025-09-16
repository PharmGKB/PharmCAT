package org.pharmgkb.pharmcat.util;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.Constants;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.phenotype.PhenotypeMap;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.PgkbGuidelineCollection;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Summary of all the genes and drugs used in PharmCAT. The expected output is a Markdown file with lists of data.
 */
public class GeneDrugSummary {
  private static final Path SUMMARY_TEMPLATE_FILE =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/util/summary.md");
  private static final Path PHENOTYPES_TEMPLATE_FILE =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/util/phenotypes.md");
  private static final String SUMMARY_FILE_NAME = "Genes-Drugs.md";
  private static final String PHENOTYPES_MD_FILE_NAME = "Phenotypes-List.md";
  private static final String PHENOTYPES_TSV_FILE_NAME = "phenotypes.tsv";
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Pattern sf_globCopyNumberPattern = Pattern.compile("\\*\\d+xN");
  private final DefinitionReader m_definitionReader;
  private final PhenotypeMap m_phenotypeMap;
  private final PgkbGuidelineCollection m_guidelineCollection;


  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("o", "output-dir", "directory to write files to", true, "directory");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      DefinitionReader definitionReader = DefinitionReader.defaultReader();

      new GeneDrugSummary(definitionReader, new PhenotypeMap(), new PgkbGuidelineCollection())
          .write(cliHelper.getValidDirectory("o", true));
    } catch (Exception e) {
      sf_logger.error("Error writing drug data", e);
    }
  }


  public GeneDrugSummary(DefinitionReader definitionReader, PhenotypeMap phenotypeMap,
      PgkbGuidelineCollection pgkbGuidelineCollection) {
    m_definitionReader = definitionReader;
    m_phenotypeMap = phenotypeMap;
    m_guidelineCollection = pgkbGuidelineCollection;
  }


  public void write(Path documentationDir) throws IOException {

    Multimap<String, DataSource> geneSourceMap = TreeMultimap.create();
    Multimap<String, PrescribingGuidanceSource> drugSourceMap = TreeMultimap.create();
    m_guidelineCollection.getGenesUsedInSource(DataSource.DPWG)
        .forEach(g -> geneSourceMap.put(g, DataSource.DPWG));
    m_guidelineCollection.getGenesUsedInSource(DataSource.CPIC)
        .forEach(g -> geneSourceMap.put(g, DataSource.CPIC));
    m_guidelineCollection.getChemicalsUsedInSource(PrescribingGuidanceSource.CPIC_GUIDELINE)
        .forEach(d -> drugSourceMap.put(d, PrescribingGuidanceSource.CPIC_GUIDELINE));
    m_guidelineCollection.getChemicalsUsedInSource(PrescribingGuidanceSource.DPWG_GUIDELINE)
        .forEach(d -> drugSourceMap.put(d, PrescribingGuidanceSource.DPWG_GUIDELINE));
    m_guidelineCollection.getChemicalsUsedInSource(PrescribingGuidanceSource.FDA_LABEL)
        .forEach(d -> drugSourceMap.put(d, PrescribingGuidanceSource.FDA_LABEL));
    m_guidelineCollection.getChemicalsUsedInSource(PrescribingGuidanceSource.FDA_ASSOC)
        .forEach(d -> drugSourceMap.put(d, PrescribingGuidanceSource.FDA_ASSOC));

    // get matcher gene list
    StringBuilder matcherGeneList = new StringBuilder()
        .append("| Gene | CPIC | DPWG |\n")
        .append("| :--- | :---: | :---: |\n");
    m_definitionReader.getGenes().stream()
        .filter(g -> !Constants.PREFER_OUTSIDE_CALL.contains(g))
        .sorted()
        .forEach(g -> appendGene(matcherGeneList, g));
    // outside call gene list
    StringBuilder outsideCallGeneList = new StringBuilder()
        .append("| Gene | CPIC | DPWG |\n")
        .append("| :--- | :---: | :---: |\n");
    Constants.PREFER_OUTSIDE_CALL.stream()
        .sorted()
        .forEach(g -> appendGene(outsideCallGeneList, g));

    // drug list
    StringBuilder drugList = new StringBuilder()
        .append("| Drug | CPIC | DPWG | FDA Label | FDA PGx Assoc |\n")
        .append("| :--- | :---: | :---: | :---: | :---: |\n");
    drugSourceMap.keySet().stream()
        .sorted()
        .forEach(d -> {
          drugList.append("| ")
              .append(d)
              .append(" | ");
          if (drugSourceMap.get(d).contains(PrescribingGuidanceSource.CPIC_GUIDELINE)) {
            drugList.append(":heavy_check_mark:");
          }
          drugList.append(" | ");
          if (drugSourceMap.get(d).contains(PrescribingGuidanceSource.DPWG_GUIDELINE)) {
            drugList.append(":heavy_check_mark:");
          }
          drugList.append(" | ");
          if (drugSourceMap.get(d).contains(PrescribingGuidanceSource.FDA_LABEL)) {
            drugList.append(":heavy_check_mark:");
          }
          drugList.append(" | ");
          if (drugSourceMap.get(d).contains(PrescribingGuidanceSource.FDA_ASSOC)) {
            drugList.append(":heavy_check_mark:");
          }
          drugList.append(" |\n");
        });

    String summaryTemplate = FileUtils.readFileToString(SUMMARY_TEMPLATE_FILE.toFile(), StandardCharsets.UTF_8);
    Path summaryFile = documentationDir.resolve(SUMMARY_FILE_NAME);
    try (BufferedWriter writer = Files.newBufferedWriter(summaryFile)) {
      writer.write(String.format(summaryTemplate, matcherGeneList, outsideCallGeneList, drugSourceMap.keySet().size(),
          drugList));
    }
    sf_logger.info("Saving summary to {}", summaryFile);


    String phenotypesTemplate = FileUtils.readFileToString(PHENOTYPES_TEMPLATE_FILE.toFile(), StandardCharsets.UTF_8);
    Path mdAlleleFile = documentationDir.resolve(PHENOTYPES_MD_FILE_NAME);
    Path tsvAlleleFile = documentationDir.resolve(PHENOTYPES_TSV_FILE_NAME);
    try (PrintWriter mdWriter = new PrintWriter(Files.newBufferedWriter(mdAlleleFile));
         PrintWriter tsvWriter = new PrintWriter(Files.newBufferedWriter(tsvAlleleFile))) {
      mdWriter.println(phenotypesTemplate);
      tsvWriter.println("Gene\tNamed Alleles\tCPIC Phenotypes\tCPIC Activity Scores\tDPWG Phenotypes");

      SortedSet<String> allGenes = new TreeSet<>(m_definitionReader.getGenes());
      m_phenotypeMap.getCpicGenes().stream().map(GenePhenotype::getGene).forEach(allGenes::add);
      m_phenotypeMap.getDpwgGenes().stream().map(GenePhenotype::getGene).forEach(allGenes::add);

      for (String gene : allGenes) {
        GenePhenotype cpicGp = m_phenotypeMap.getPhenotype(gene, DataSource.CPIC);
        GenePhenotype dpwgGp = m_phenotypeMap.getPhenotype(gene, DataSource.DPWG);

        SortedSet<String> haplotypes = new TreeSet<>(new HaplotypeNameComparator());
        if (Constants.PREFER_OUTSIDE_CALL.contains(gene)) {
          if (cpicGp != null) {
            haplotypes.addAll(cpicGp.getHaplotypes().keySet());
          }
          if (dpwgGp != null) {
            haplotypes.addAll(dpwgGp.getHaplotypes().keySet());
          }

        } else {
          m_definitionReader.getDefinitionFile(gene).getNamedAlleles().stream()
              .map(NamedAllele::getName)
              .forEach(haplotypes::add);
        }

        SortedSet<String> cpicPhenotypes = new TreeSet<>();
        SortedSet<String> cpicScores = new TreeSet<>();
        collectDiplotypeMetadata(cpicGp, cpicPhenotypes, cpicScores);
        SortedSet<String> dpwgPhenotypes = new TreeSet<>();
        SortedSet<String> dpwgScores = new TreeSet<>();
        collectDiplotypeMetadata(dpwgGp, dpwgPhenotypes, dpwgScores);

        // we use semicolon as a separator so make sure it's not used in values
        if (haplotypes.stream().anyMatch(v -> v.contains(";"))) {
          throw new IllegalStateException(gene + " haplotypes has semicolon: " +
              haplotypes.stream()
                  .filter(v -> v.contains(";"))
                  .map(v -> "\"" + v + "\"")
                  .collect(Collectors.joining(", ")));
        }
        if (cpicPhenotypes.stream().anyMatch(v -> v.contains(";"))) {
          throw new IllegalStateException(gene + " CPIC phenotypes has semicolon");
        }
        if (cpicPhenotypes.stream().anyMatch(v -> v.contains("metabolizer"))) {
          throw new IllegalStateException(gene + " CPIC phenotypes has lower cased \"Metabolizer\"");
        }
        if (cpicPhenotypes.stream().anyMatch(v -> v.toLowerCase().contains("metabolise"))) {
          throw new IllegalStateException(gene + " CPIC phenotype uses metaboliSe");
        }
        if (cpicScores.stream().anyMatch(v -> v.contains(";"))) {
          throw new IllegalStateException(gene + " CPIC scores has semicolon");
        }
        if (dpwgPhenotypes.stream().anyMatch(v -> v.contains(";"))) {
          throw new IllegalStateException(gene + " DPWG phenotypes has semicolon");
        }
        if (dpwgPhenotypes.stream().anyMatch(v -> v.toLowerCase().contains("metabolise"))) {
          throw new IllegalStateException(gene + " DPWG phenotype uses metaboliSe");
        }

        if (!cpicPhenotypes.isEmpty() && !dpwgPhenotypes.isEmpty()) {
          if (!(cpicPhenotypes.containsAll(dpwgPhenotypes) || dpwgPhenotypes.containsAll(cpicPhenotypes))) {
            System.out.println("WARNING: " + gene + " has CPIC/DPWG phenotype mismatch");
            System.out.println("         CPIC: " + cpicPhenotypes);
            System.out.println("         DPWG: " + dpwgPhenotypes);
          }
        }

        String type = Constants.isVariantGene(gene) ? "Named Variants" : "Named Alleles";
        mdWriter.println("### " + gene);
        if (!m_guidelineCollection.getGenesWithRecommendations().contains(gene)) {
          mdWriter.println("<p>No recommendations available for this gene.</p>");
        }
        mdWriter.println("<table>");
        mdWriter.println("<tr>");
        mdWriter.println("<th style=\"text-align: left\">" + type + "</th>");
        if (!cpicPhenotypes.isEmpty()) {
          mdWriter.println("<th style=\"text-align: left\">CPIC Phenotypes</th>");
        }
        if (!cpicScores.isEmpty()) {
          mdWriter.println("<th style=\"text-align: left\">CPIC Activity Scores</th>");
        }
        if (!dpwgPhenotypes.isEmpty()) {
          mdWriter.println("<th style=\"text-align: left\">DPWG Phenotypes</th>");
        }
        if (!dpwgScores.isEmpty()) {
          mdWriter.println("<th style=\"text-align: left\">DPWG Activity Scores</th>");
        }
        mdWriter.println("</tr>");
        mdWriter.println("<tr>");
        mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
        haplotypes.stream()
            .filter(v -> !sf_globCopyNumberPattern.matcher(v).matches())
            .forEach(v -> mdWriter.println("<li>" + v + "</li>"));
        mdWriter.println("</ul></td>");

        if (!cpicPhenotypes.isEmpty()) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          cpicPhenotypes.forEach(v -> mdWriter.println("<li>" + v + "</li>"));
          mdWriter.println("</ul></td>");
        }
        if (!cpicScores.isEmpty()) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          cpicScores.forEach(v -> mdWriter.println("<li>" + v + "</li>"));
          mdWriter.println("</ul></td>");
        }
        if (!dpwgPhenotypes.isEmpty()) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          dpwgPhenotypes.forEach(v -> mdWriter.println("<li>" + v + "</li>"));
          mdWriter.println("</ul></td>");
        }
        if (!dpwgScores.isEmpty()) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          dpwgScores.forEach(v -> mdWriter.println("<li>" + v + "</li>"));
          mdWriter.println("</ul></td>");
        }
        mdWriter.println("</tr>");
        mdWriter.println("</table>");
        mdWriter.println();

        tsvWriter.print(gene);
        tsvWriter.print("\t");
        tsvWriter.print(String.join(";", haplotypes));
        tsvWriter.print("\t");
        tsvWriter.print(String.join(";", cpicPhenotypes));
        tsvWriter.print("\t");
        tsvWriter.print(String.join(";", cpicScores));
        tsvWriter.print("\t");
        tsvWriter.print(String.join(";", dpwgPhenotypes));
        tsvWriter.println();
      }


      mdWriter.println("## Allele Presence Genes");
      mdWriter.println();
      mdWriter.println("""
      PharmCAT provides recommendations for these genes based on the presence or absence of specific named alleles.
      """);
      mdWriter.println();

      for (String gene : Constants.ALLELE_PRESENCE_GENES) {
        mdWriter.println("### " + gene);

        if (gene.equals("HLA-A")) {
          Constants.HLA_A_ALLELES.forEach(a -> mdWriter.println("* " + a));

        } else if (gene.equals("HLA-B")) {
          Constants.HLA_B_ALLELES.forEach(a -> mdWriter.println("* " + a));
        } else {
          throw new IllegalStateException("Don't know how to handle " + gene);
        }
        mdWriter.println();

        tsvWriter.print(gene);
        tsvWriter.print("\t");
        if (gene.equals("HLA-A")) {
          tsvWriter.print(String.join(";", Constants.HLA_A_ALLELES));
        } else {
          tsvWriter.print(String.join(";", Constants.HLA_B_ALLELES));
        }
        tsvWriter.println("\t\t\t");
      }
    }
    sf_logger.info("Saving allele data to {}", summaryFile);


    Path nat2ExceptionFile = documentationDir.resolve("methods/unphasedPriorities-NAT2.md");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(nat2ExceptionFile))) {
      DefinitionExemption exemption = Objects.requireNonNull(m_definitionReader.getExemption("NAT2"));
      AtomicInteger count = new AtomicInteger(0);
      int half = (exemption.getUnphasedDiplotypePriorities().size() + 1) / 2;
      List<String> col1 = new ArrayList<>();
      List<String> col2 = new ArrayList<>();
      exemption.getUnphasedDiplotypePriorities().forEach(entry -> {
        StringBuilder builder = new StringBuilder();
        builder.append("* __")
            .append(entry.getPick().replaceAll("\\*", "\\\\*"))
            .append("__");
        for (String dip : entry.getList()) {
          if (!dip.equals(entry.getPick())) {
            builder.append(", ")
                .append(dip.replaceAll("\\*", "\\\\*"));
          }
        }
        List<String> col = count.getAndIncrement() > half ? col2 : col1;
        col.add(builder.toString());
      });
      for (int x = 0; x < col1.size(); x++) {
        writer.print("| " + col1.get(x) + " | ");
        if (col2.size() > x) {
          writer.print(col2.get(x));
        }
        writer.println(" |");
      }
    }
  }

  private void collectDiplotypeMetadata(GenePhenotype dpwgGp, SortedSet<String> dpwgPhenotypes, SortedSet<String> dpwgScores) {
    if (dpwgGp != null) {
      dpwgGp.getDiplotypes().forEach(d -> {
        if (d.getGeneResult() != null) {
          dpwgPhenotypes.add(d.getGeneResult());
        }
        String score = d.getActivityScore();
        if (score != null && !score.equals(TextConstants.NA)) {
          dpwgScores.add(d.getActivityScore());
        }
      });
    }
  }

  private void printPhenotype(PrintWriter writer, GenePhenotype gp, String haplotype) {
    if (gp != null) {
      writer.print(StringUtils.stripToEmpty(gp.getHaplotypeFunction(haplotype)));
      writer.print("\t");
      writer.print(StringUtils.stripToEmpty(gp.getHaplotypeActivity(haplotype)));
    } else {
      writer.print("\t");
    }
    writer.print("\t");
  }


  private void appendGene(StringBuilder builder, String gene) {
    builder.append("| [")
        .append(gene)
        .append("](/Phenotypes-List#")
        .append(gene.toLowerCase())
        .append(") | ");
    if (m_guidelineCollection.getGenesUsedInSource(DataSource.CPIC).contains(gene)) {
      builder.append(":heavy_check_mark:");
    }
    builder.append(" | ");
    if (m_guidelineCollection.getGenesUsedInSource(DataSource.DPWG).contains(gene)) {
      builder.append(":heavy_check_mark:");
    }
    builder.append(" |\n");
  }
}
