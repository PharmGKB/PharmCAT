package org.pharmgkb.pharmcat.util;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import com.google.common.base.Charsets;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.phenotype.PhenotypeMap;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.PgkbGuidelineCollection;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Summary of all the genes and drugs used in PharmCAT. Expected output is a markdown file with lists of data.
 */
public class GeneDrugSummary {
  private static final Path SUMMARY_TEMPLATE_FILE =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/util/summary.md");
  private static final Path PHENOTYPES_TEMPLATE_FILE =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/util/phenotypes.md");
  private static final String SUMMARY_FILE_NAME = "Genes-and-Drugs.md";
  private static final String PHENOTYPES_MD_FILE_NAME = "Phenotypes-List.md";
  private static final String PHENOTYPES_TSV_FILE_NAME = "phenotypes.tsv";
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Set<String> PREFER_OUTSIDE_CALL = ImmutableSet.of("CYP2D6", "F5", "HLA-A", "HLA-B", "MT-RNR1");
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
    Multimap<String, DataSource> drugSourceMap = TreeMultimap.create();
    m_guidelineCollection.getGenesUsedInSource(DataSource.DPWG)
        .forEach(g -> geneSourceMap.put(g, DataSource.DPWG));
    m_guidelineCollection.getGenesUsedInSource(DataSource.CPIC)
        .forEach(g -> geneSourceMap.put(g, DataSource.CPIC));
    m_guidelineCollection.getChemicalsUsedInSource(DataSource.CPIC)
        .forEach(d -> drugSourceMap.put(d, DataSource.CPIC));
    m_guidelineCollection.getChemicalsUsedInSource(DataSource.DPWG)
        .forEach(d -> drugSourceMap.put(d, DataSource.DPWG));

    // get matcher gene list
    StringBuilder matcherGeneList = new StringBuilder()
        .append("| Gene | CPIC | PharmGKB-DPWG |\n")
        .append("| :--- | :---: | :---: |\n");
    m_definitionReader.getGenes().stream()
        .filter(g -> !PREFER_OUTSIDE_CALL.contains(g))
        .sorted()
        .forEach(g -> appendGene(matcherGeneList, g));
    // outside call gene list
    StringBuilder outsideCallGeneList = new StringBuilder()
        .append("| Gene | CPIC | PharmGKB-DPWG |\n")
        .append("| :--- | :---: | :---: |\n");
    PREFER_OUTSIDE_CALL.stream()
        .sorted()
        .forEach(g -> appendGene(outsideCallGeneList, g));

    // drug list
    StringBuilder drugList = new StringBuilder()
        .append("| Drug | CPIC | PharmGKB-DPWG |\n")
        .append("| :--- | :---: | :---: |\n");
    drugSourceMap.keySet().stream()
        .sorted()
        .forEach(d -> {
          drugList.append("| ")
              .append(d)
              .append(" | ");
          if (drugSourceMap.get(d).contains(DataSource.CPIC)) {
            drugList.append(":heavy_check_mark:");
          }
          drugList.append(" | ");
          if (drugSourceMap.get(d).contains(DataSource.DPWG)) {
            drugList.append(":heavy_check_mark:");
          }
          drugList.append(" |\n");
        });

    String summaryTemplate = FileUtils.readFileToString(SUMMARY_TEMPLATE_FILE.toFile(), Charsets.UTF_8);
    String phenotypesTemplate = FileUtils.readFileToString(PHENOTYPES_TEMPLATE_FILE.toFile(), Charsets.UTF_8);
    Path summaryFile = documentationDir.resolve(SUMMARY_FILE_NAME);
    try (BufferedWriter writer = Files.newBufferedWriter(summaryFile)) {
      writer.write(String.format(summaryTemplate, matcherGeneList, outsideCallGeneList, drugList));
    }
    sf_logger.info("Saving summary to {}", summaryFile);

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
        if (PREFER_OUTSIDE_CALL.contains(gene)) {
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

        // we use semicolon as separator so make sure it's not used in values
        if (haplotypes.stream().anyMatch(v -> v.contains(";"))) {
          throw new IllegalStateException(gene + " haplotypes has comma");
        }
        if (cpicPhenotypes.stream().anyMatch(v -> v.contains(";"))) {
          throw new IllegalStateException(gene + " CPIC phenotypes has comma");
        }
        if (cpicPhenotypes.stream().anyMatch(v -> v.contains("metabolizer"))) {
          throw new IllegalStateException(gene + " CPIC phenotypes has lower cased \"Metabolizer\"");
        }
        if (cpicPhenotypes.stream().anyMatch(v -> v.toLowerCase().contains("metabolise"))) {
          throw new IllegalStateException(gene + " CPIC phenotype uses metaboliSe");
        }
        if (cpicScores.stream().anyMatch(v -> v.contains(";"))) {
          throw new IllegalStateException(gene + " CPIC scores has comma");
        }
        if (dpwgPhenotypes.stream().anyMatch(v -> v.contains(";"))) {
          throw new IllegalStateException(gene + " DPWG phenotypes has comma");
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

        boolean isNamedAllele = haplotypes.first().startsWith("*") && !gene.equals("UGT1A1");
        String type = isNamedAllele ? "Named Alleles" : "Variants";
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
        mdWriter.print(haplotypes.stream()
            .filter(v -> !sf_globCopyNumberPattern.matcher(v).matches())
            .map(v -> "<li>" + v + "</li>").collect(Collectors.joining()));
        mdWriter.println("</ul></td>");

        if (!cpicPhenotypes.isEmpty()) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          mdWriter.print(cpicPhenotypes.stream().map(v -> "<li>" + v + "</li>").collect(Collectors.joining()));
          mdWriter.println("</ul></td>");
        }
        if (!cpicScores.isEmpty()) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          mdWriter.print(cpicScores.stream().map(v -> "<li>" + v + "</li>").collect(Collectors.joining()));
          mdWriter.println("</ul></td>");
        }
        if (!dpwgPhenotypes.isEmpty()) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          mdWriter.print(dpwgPhenotypes.stream().map(v -> "<li>" + v + "</li>").collect(Collectors.joining()));
          mdWriter.println("</ul></td>");
        }
        if (!dpwgScores.isEmpty()) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          mdWriter.print(dpwgScores.stream().map(v -> "<li>" + v + "</li>").collect(Collectors.joining()));
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
    }
    sf_logger.info("Saving allele data to {}", summaryFile);
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
