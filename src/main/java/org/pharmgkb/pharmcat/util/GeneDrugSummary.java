package org.pharmgkb.pharmcat.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.lang.invoke.MethodHandles;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.io.IOUtils;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.reporter.DrugCollection;
import org.pharmgkb.pharmcat.reporter.PgkbGuidelineCollection;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Summary of all the genes and drugs used in PharmCAT. Expected output is a markdown file with lists of data.
 */
public class GeneDrugSummary {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final String SUMMARY_REPORT = "Genes-and-Drugs.md";
  private static final String TEMPLATE_FILE_NAME = "summary.md";
  private static final Set<String> PREFER_OUTSIDE_CALL = ImmutableSet.of("CYP2D6", "G6PD", "HLA-A", "HLA-B", "MT-RNR1");

  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("o", "output-dir", "directory to write files to", true, "directory");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      DefinitionReader definitionReader = new DefinitionReader();
      definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);

      new GeneDrugSummary().write(
          cliHelper.getValidDirectory("o", true),
          definitionReader.getGeneAlleleCount(),
          new DrugCollection(),
          new PhenotypeMap()
      );
    } catch (Exception e) {
      sf_logger.error("Error writing drug data", e);
    }
  }

  public void write(Path documentationDir, Map<String,Integer> geneAlleleCount, DrugCollection drugs, PhenotypeMap phenotypeMap) throws IOException {
    PgkbGuidelineCollection pgkbGuidelineCollection = new PgkbGuidelineCollection();
    Multimap<String, DataSource> geneToUsageMap = makeGeneToUsageMap(drugs, pgkbGuidelineCollection);
    try (
        InputStream mdStream = getClass().getResourceAsStream(TEMPLATE_FILE_NAME);
        StringWriter templateWriter = new StringWriter()
    ) {
      if (mdStream == null) {
        throw new RuntimeException("No summary report template found");
      }
      IOUtils.copy(mdStream, templateWriter, Charset.defaultCharset());
      String mdTemplate = templateWriter.toString();

      File summaryFile = documentationDir.resolve(SUMMARY_REPORT).toFile();
      try (FileWriter fw = new FileWriter(summaryFile)) {
        String matcherGeneList = geneAlleleCount.keySet().stream()
            .filter(g -> !PREFER_OUTSIDE_CALL.contains(g))
            .sorted().map(g -> "- " + g + " (" + geneAlleleCount.get(g) + " alleles)" + makeUsageSuffix(geneToUsageMap.get(g)))
            .collect(Collectors.joining("\n"));

        String outsideGeneList = PREFER_OUTSIDE_CALL.stream()
            .map(g ->
              phenotypeMap.lookup(g).map((genePhenotype) -> "- " + genePhenotype.getGene() + " (" + genePhenotype.getHaplotypes().size() + " alleles)" + makeUsageSuffix(geneToUsageMap.get(g)))
                  .orElse("- " + g + makeUsageSuffix(geneToUsageMap.get(g))
              )
            )
            .collect(Collectors.joining("\n"));

        String drugList = String.join("\n", makeDrugListItems(drugs, pgkbGuidelineCollection));

        IOUtils.write(String.format(mdTemplate, matcherGeneList, outsideGeneList, drugList), fw);
      }
      sf_logger.info("Saving summary file to {}", summaryFile);
    }
  }

  private Multimap<String,DataSource> makeGeneToUsageMap(DrugCollection drugs, PgkbGuidelineCollection pgkbGuidelineCollection) {
    Multimap<String,DataSource> geneMap = TreeMultimap.create();
    drugs.forEach((drug) -> {
      for (String gene : drug.getGenes()) {
        geneMap.put(gene, drug.getSource());
      }
    });
    pgkbGuidelineCollection.getGenes().forEach(g -> geneMap.put(g, DataSource.DPWG));
    return geneMap;
  }

  private List<String> makeDrugListItems(DrugCollection drugCollection, PgkbGuidelineCollection pgkbGuidelineCollection) {
    Multimap<String,DataSource> drugSourceMap = TreeMultimap.create();
    for (Drug drug : drugCollection) {
      drugSourceMap.put(drug.getDrugName(), drug.getSource());
    }

    for (GuidelinePackage guidelinePackage : pgkbGuidelineCollection.getGuidelinePackages()) {
      for (String drugName : guidelinePackage.getDrugs()) {
        drugSourceMap.put(drugName, DataSource.DPWG);
      }
    }

    return drugSourceMap.keySet().stream()
        .map(name -> "- " + name + makeUsageSuffix(drugSourceMap.get(name)))
        .collect(Collectors.toList());
  }

  private String makeUsageSuffix(Collection<DataSource> sources) {
    if (sources == null || sources.size() == 0) {
      return "";
    } else {
      return " (" + sources.stream().map(DataSource::getDisplayName).collect(Collectors.joining(", ")) + ")";
    }
  }
}
