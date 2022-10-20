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
import org.pharmgkb.pharmcat.reporter.DrugCollection;
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
  // TODO: revert when HLA's are supported again
  //private static final Set<String> PREFER_OUTSIDE_CALL = ImmutableSet.of("CYP2D6", "F5", "HLA-A", "HLA-B", "MT-RNR1");
  private static final Set<String> PREFER_OUTSIDE_CALL = ImmutableSet.of("CYP2D6", "F5", "MT-RNR1");
  private final DefinitionReader m_definitionReader;
  private final PhenotypeMap m_phenotypeMap;
  private final DrugCollection m_cpicDrugs;
  private final PgkbGuidelineCollection m_dpwgDrugs;


  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("o", "output-dir", "directory to write files to", true, "directory");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      DefinitionReader definitionReader = new DefinitionReader();
      definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);

      new GeneDrugSummary(definitionReader, new PhenotypeMap(), new DrugCollection(), new PgkbGuidelineCollection())
          .write(cliHelper.getValidDirectory("o", true));
    } catch (Exception e) {
      sf_logger.error("Error writing drug data", e);
    }
  }


  public GeneDrugSummary(DefinitionReader definitionReader, PhenotypeMap phenotypeMap, DrugCollection drugCollection,
      PgkbGuidelineCollection pgkbGuidelineCollection) {
    m_definitionReader = definitionReader;
    m_phenotypeMap = phenotypeMap;
    m_cpicDrugs = drugCollection;
    m_dpwgDrugs = pgkbGuidelineCollection;
  }


  public void write(Path documentationDir) throws IOException {

    Multimap<String, DataSource> geneSourceMap = TreeMultimap.create();
    Multimap<String, DataSource> drugSourceMap = TreeMultimap.create();
    m_cpicDrugs.forEach((drug) -> {
      drugSourceMap.put(drug.getDrugName(), DataSource.CPIC);
      for (String gene : drug.getGenes()) {
        geneSourceMap.put(gene, DataSource.CPIC);
      }
    });
    m_dpwgDrugs.getGenes()
        .forEach(g -> geneSourceMap.put(g, DataSource.DPWG));
    m_dpwgDrugs.getChemicals()
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

    // get drug list
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
      tsvWriter.println("Gene\tNamed Alleles\tCPIC Phenotypes\tCPIC Activity Scores\tDPWG Phenotyeps");

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
        SortedSet<String> dpwgPhenotypes = new TreeSet<>();
        if (cpicGp != null) {
          cpicGp.getDiplotypes().forEach(d -> {
            if (d.getGeneResult() != null) {
              cpicPhenotypes.add(d.getGeneResult());
            }
            String score = d.getLookupKey();
            if (score != null && !score.equals(TextConstants.NA) && !score.equals(d.getGeneResult())) {
              cpicScores.add(d.getLookupKey());
            }
          });
        }
        if (dpwgGp != null) {
          dpwgGp.getDiplotypes().forEach(d -> {
            if (d.getGeneResult() != null) {
              dpwgPhenotypes.add(d.getGeneResult());
            }
          });
        }

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

        if (cpicPhenotypes.size() > 0 && dpwgPhenotypes.size() > 0) {
          if (!(cpicPhenotypes.containsAll(dpwgPhenotypes) || dpwgPhenotypes.containsAll(cpicPhenotypes))) {
            System.out.println("WARNING: " + gene + " has CPIC/DPWG phenotype mismatch");
            System.out.println("         CPIC: " + cpicPhenotypes);
            System.out.println("         DPWG: " + dpwgPhenotypes);
          }
        }

        boolean isNamedAllele = haplotypes.first().startsWith("*") && !gene.equals("UGT1A1");
        String type = isNamedAllele ? "Named Alleles" : "Variants";
        mdWriter.println("### " + gene);
        mdWriter.println("<table>");
        mdWriter.println("<tr>");
        mdWriter.println("<th style=\"text-align: left\">" + type + "</th>");
        if (cpicPhenotypes.size() > 0) {
          mdWriter.println("<th style=\"text-align: left\">CPIC Phenotypes</th>");
        }
        if (cpicScores.size() > 0) {
          mdWriter.println("<th style=\"text-align: left\">CPIC Activity Scores</th>");
        }
        if (dpwgPhenotypes.size() > 0) {
          mdWriter.println("<th style=\"text-align: left\">DPWG Phenotypes</th>");
        }
        mdWriter.println("</tr>");
        mdWriter.println("<tr>");
        mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
        mdWriter.print(haplotypes.stream().map(v -> "<li>" + v + "</li>").collect(Collectors.joining()));
        mdWriter.println("</ul></td>");

        if (cpicPhenotypes.size() > 0) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          mdWriter.print(cpicPhenotypes.stream().map(v -> "<li>" + v + "</li>").collect(Collectors.joining()));
          mdWriter.println("</ul></td>");
        }
        if (cpicScores.size() > 0) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          mdWriter.print(cpicScores.stream().map(v -> "<li>" + v + "</li>").collect(Collectors.joining()));
          mdWriter.println("</ul></td>");
        }
        if (dpwgPhenotypes.size() > 0) {
          mdWriter.print("<td style=\"vertical-align: top\"><ul style=\"padding-left: 1rem\">");
          mdWriter.print(dpwgPhenotypes.stream().map(v -> "<li>" + v + "</li>").collect(Collectors.joining()));
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
    if (m_phenotypeMap.getPhenotype(gene, DataSource.CPIC) != null ||
        m_cpicDrugs.getAllReportableGenes().contains(gene)) {
      builder.append(":heavy_check_mark:");
    }
    builder.append(" | ");
    if (m_phenotypeMap.getPhenotype(gene, DataSource.DPWG) != null ||
        m_dpwgDrugs.getGenes().contains(gene)) {
      builder.append(":heavy_check_mark:");
    }
    builder.append(" |\n");
  }
}
