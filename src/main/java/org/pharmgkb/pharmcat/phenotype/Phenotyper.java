package org.pharmgkb.pharmcat.phenotype;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.definition.MessageList;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.definition.ReferenceAlleleMap;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.Reporter;
import org.pharmgkb.pharmcat.reporter.io.OutsideCallParser;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.DataManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class takes genotyping information from the {@link NamedAlleleMatcher} and from outside allele calls then
 * interprets them into phenotype assignments, diplotype calls, and other information needed for subsequent use in the
 * {@link Reporter}. The data is compiled into {@link GeneReport} objects which can then serialized
 */
public class Phenotyper {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private final SortedSet<GeneReport> f_geneReports = new TreeSet<>();

  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        //
        .addOption("vcf", "sample-vcf", "input sample file (VCF)", false, "vcf")
        .addOption("c", "call-file", "named allele call JSON file", false, "call-file-path")
        // optional-data
        .addOption("a", "outside-call-file", "optional, outside call TSV file", false, "outside-file-path")
        // output
        .addOption("f", "output-file", "file path to write Phenotyper JSON data to", true, "output-file-path")
        // research
        .addOption("r", "research-mode", "enable research mode")
        .addOption("cyp2d6", "research-cyp2d6", "call CYP2D6 (must also use research mode)")
        ;
    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }
      Path vcfFile = cliHelper.hasOption("vcf") ? cliHelper.getValidFile("vcf", true) : null;
      Path callFile = cliHelper.hasOption("c") ? cliHelper.getValidFile("c", true) : null;
      Path outsideCallPath = cliHelper.hasOption("a") ? cliHelper.getValidFile("a", true) : null;
      Path outputFile = cliHelper.getPath("f");

      Preconditions.checkArgument(callFile != null ^ vcfFile != null, "Can use VCF file or Matcher JSON file, not both");
      Preconditions.checkArgument(callFile == null || Files.isRegularFile(callFile), "Call file does not exist or is not a regular file");
      Preconditions.checkArgument(vcfFile == null || Files.isRegularFile(vcfFile), "Sample VCF file does not exist or is not a regular file");

      List<GeneCall> calls;
      Map<String,Collection<String>> variantWarnings = new HashMap<>();
      if (callFile != null) {
        calls = new ResultSerializer().fromJson(callFile).getGeneCalls();
      } else {
        DefinitionReader definitionReader = new DefinitionReader();
        definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
        boolean callCyp2d6 = cliHelper.hasOption("r") && cliHelper.hasOption("cyp2d6");
        NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(definitionReader, false, false, callCyp2d6);
        Result result = namedAlleleMatcher.call(vcfFile);
        calls = result.getGeneCalls();
        variantWarnings = result.getVcfWarnings();
      }

      //Load the outside calls if it's available
      List<OutsideCall> outsideCalls = new ArrayList<>();
      if (outsideCallPath != null) {
        Preconditions.checkArgument(Files.exists(outsideCallPath));
        Preconditions.checkArgument(Files.isRegularFile(outsideCallPath));
        outsideCalls = OutsideCallParser.parse(outsideCallPath);
      }

      Phenotyper phenotyper = new Phenotyper(calls, outsideCalls, variantWarnings);
      phenotyper.write(outputFile);
    } catch (IOException e) {
      e.printStackTrace();
      System.exit(1);
    }
  }

  /**
   * Public constructor. This needs {@link GeneCall} objects from the {@link NamedAlleleMatcher} and {@link OutsideCall}
   * objects coming from other allele calling sources. This relies on reading definition files as well.
   * @param geneCalls a List of {@link GeneCall} objects
   * @param outsideCalls a List of {@link OutsideCall} objects
   */
  public Phenotyper(List<GeneCall> geneCalls, List<OutsideCall> outsideCalls, @Nullable Map<String, Collection<String>> variantWarnings) {
    ReferenceAlleleMap referenceAlleleMap = new ReferenceAlleleMap();
    PhenotypeMap phenotypeMap = new PhenotypeMap();

    for (GeneCall geneCall : geneCalls) {
      GeneReport geneReport = new GeneReport(geneCall);
      DiplotypeFactory diplotypeFactory = new DiplotypeFactory(
          geneReport.getGene(),
          phenotypeMap.lookup(geneReport.getGene()).orElse(null),
          referenceAlleleMap.get(geneReport.getGene()));
      geneReport.setDiplotypes(diplotypeFactory, geneCall);

      f_geneReports.add(geneReport);
    }

    for (OutsideCall outsideCall : outsideCalls) {
      GeneReport geneReport = findGeneReport(outsideCall.getGene())
          .orElse(null);

      if (geneReport != null) {
        // DO NOT allow both a report from the matcher/caller and from an outside source
        if (geneReport.getCallSource() == CallSource.MATCHER && geneReport.isCalled()) {
          throw new ParseException("Cannot specify outside call for " + geneReport.getGene() + ", it is already called in sample data");
        } else {
          // if the caller already made a report but it's uncalled let's remove it so we can replace it
          removeGeneReport(outsideCall.getGene());
        }
      }
      geneReport = new GeneReport(outsideCall);
      f_geneReports.add(geneReport);

      DiplotypeFactory diplotypeFactory = new DiplotypeFactory(
          geneReport.getGene(),
          phenotypeMap.lookup(geneReport.getGene()).orElse(null),
          referenceAlleleMap.get(geneReport.getGene()));
      geneReport.setDiplotypes(diplotypeFactory, outsideCall);

      f_geneReports.add(geneReport);
    }

    f_geneReports.forEach(r -> r.addVariantWarningMessages(variantWarnings));

    try {
      MessageList messageList = new MessageList();
      getGeneReports().forEach(messageList::addMatchingMessagesTo);
    } catch (IOException ex) {
      throw new RuntimeException("Could not load messages", ex);
    }
  }

  /**
   * Write the collection of {@link GeneReport} objects
   * @param outputPath the path to write a JSON file of data to
   * @throws IOException can occur from writing the file to the filesystem
   */
  public void write(Path outputPath) throws IOException {
    try (BufferedWriter writer = Files.newBufferedWriter(outputPath, StandardCharsets.UTF_8)) {
      Gson gson = new GsonBuilder().serializeNulls().excludeFieldsWithoutExposeAnnotation()
          .setPrettyPrinting().create();
      writer.write(gson.toJson(f_geneReports));
      sf_logger.info("Writing Phenotyper JSON to " + outputPath);
    }
  }

  /**
   * Get the {@link GeneReport} objects that were created from call data in the constructor
   * @return a SortedSet of {@link GeneReport} objects
   */
  public SortedSet<GeneReport> getGeneReports() {
    return f_geneReports;
  }

  /**
   * Find a {@link GeneReport} based on the gene symbol
   * @param geneSymbol a gene symbol
   */
  public Optional<GeneReport> findGeneReport(String geneSymbol) {
    return getGeneReports().stream().filter(r -> r.getGene().equals(geneSymbol)).findFirst();
  }

  private void removeGeneReport(String geneSymbol) {
    findGeneReport(geneSymbol).ifPresent(f_geneReports::remove);
  }

  /**
   * Read gene resport information from a given JSON file path. This should be the JSON output of this class.
   * @param filePath a path to an existing JSON file
   * @return a List of {@link GeneReport} objects
   * @throws IOException can occur if file is unable to be read
   */
  public static List<GeneReport> readGeneReports(Path filePath) throws IOException {
    Preconditions.checkNotNull(filePath);
    Preconditions.checkArgument(filePath.toFile().exists());
    Preconditions.checkArgument(filePath.toFile().isFile());
    Gson gson = new GsonBuilder().create();
    try (BufferedReader reader = Files.newBufferedReader(filePath)) {
      return Arrays.asList(gson.fromJson(reader, GeneReport[].class));
    }
  }
}
