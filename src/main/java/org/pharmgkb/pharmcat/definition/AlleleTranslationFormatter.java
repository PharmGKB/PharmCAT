package org.pharmgkb.pharmcat.definition;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.stream.Collectors;
import org.pharmgkb.pharmcat.definition.model.AlleleTranslation;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class will take allele definition files in the <code>translations</code> directory, reformat them, and write
 * them out to the <code>output</code> directory. The new format is more readily parsed for the PharmCAT haplotyper.
 *
 * @author Ryan Whaley
 * @author Alex Frase
 */
public class AlleleTranslationFormatter {
  private static final Logger sf_logger = LoggerFactory.getLogger(AlleleTranslationFormatter.class);
  private static final SimpleDateFormat sf_dateFormat = new SimpleDateFormat("MM/dd/yy");
  private static final String OUTPUT_FORMAT_VERSION = "1";
  private static final String sf_outputDir = "output";

  public static void main(String[] args) {
    AlleleTranslationFormatter formatter = new AlleleTranslationFormatter();

    try {
      formatter.run();
    }
    catch (IOException ex) {
      sf_logger.error("Couldn't run formatter", ex);
      System.exit(1);
    }

    System.exit(0);
  }

  /**
   * Only grab .tsv files
   */
  private static DirectoryStream.Filter<Path> onlyTsvs = p -> p.toString().endsWith(".tsv");

  /**
   * Does the work for stepping through the files and applying the format
   * @throws IOException can occur due to Filesystem errors
   */
  private void run() throws IOException {
    Path translationsDirectory = Paths.get("translations");

    if (!translationsDirectory.toFile().exists() || !translationsDirectory.toFile().isDirectory()) {
      throw new IOException(translationsDirectory + " is not a directory");
    }

    Files.newDirectoryStream(translationsDirectory, onlyTsvs).forEach(f -> {
      sf_logger.info("Parsing {}", f);

      AlleleTranslation alleleTranslation = new AlleleTranslationParser(f).parse();
      try (BufferedWriter output = new BufferedWriter(
          new OutputStreamWriter(Files.newOutputStream(
              Paths.get(sf_outputDir, alleleTranslation.getGeneSymbol() + ".tsv")), StandardCharsets.UTF_8))) {

        writeHeader(output, alleleTranslation);
        writeAlleles(output, alleleTranslation);
        writeNotes(output, alleleTranslation);

      } catch (IOException e) {
        throw new RuntimeException("Couldn't write output", e);
      }

      sf_logger.info("Parsed sheet for {}", alleleTranslation.getGeneSymbol());
    });
  }

  /**
   * Write the new header section to the supplied output
   * @param output the stream to write to
   * @param alleleTranslation the allele translation to read from
   * @throws IOException
   */
  private void writeHeader(BufferedWriter output, AlleleTranslation alleleTranslation) throws IOException {
    output.write("FormatVersion\t"    + OUTPUT_FORMAT_VERSION + "\n");
    output.write("GeneName\t"         + alleleTranslation.getGeneSymbol() + "\n");
    output.write("GeneRefSeq\t"       + alleleTranslation.getRefSeqGene() + "\n");
    output.write("GeneOrientation\t"  + alleleTranslation.getOrientation() + "\n");
    output.write("ContentDate\t"      + sf_dateFormat.format(alleleTranslation.getModificationDate()) + "\n");
    output.write("ContentVersion\t"   + alleleTranslation.getContentVersion() + "\n");
    output.write("GenomeBuild\t"      + alleleTranslation.getGenomeBuild() + "\n");
    output.write("ChrName\t"          + alleleTranslation.getChromoName() + "\n");
    output.write("ChrRefSeq\t"        + alleleTranslation.getRefSeqChromo() + "\n");
    output.write("ProteinRefSeq\t"    + alleleTranslation.getRefSeqProtein() + "\n");
    output.write("NumVariants\t"      + alleleTranslation.getVariants().length + "\n");

    output.write("ResourceNote\t\t\t\t");
    output.write(Arrays.stream(alleleTranslation.getVariants()).map(VariantLocus::getResourceNote).collect(Collectors.joining("\t")));
    output.write("\n");

    output.write("ProteinNote\t\t\t\t");
    output.write(Arrays.stream(alleleTranslation.getVariants()).map(VariantLocus::getProteinNote).collect(Collectors.joining("\t")));
    output.write("\n");

    output.write("ChrPosition\t\t\t\t");
    output.write(Arrays.stream(alleleTranslation.getVariants()).map(VariantLocus::getChrPosition).collect(Collectors.joining("\t")));
    output.write("\n");

    output.write("GenePosition\t\t\t\t");
    output.write(Arrays.stream(alleleTranslation.getVariants()).map(VariantLocus::getGenePosition).collect(Collectors.joining("\t")));
    output.write("\n");

    output.write("rsID\t\t\t\t");
    output.write(Arrays.stream(alleleTranslation.getVariants()).map(VariantLocus::getRsid).collect(Collectors.joining("\t")));
    output.write("\n");

    output.write("Header\tID\tName\tFunctionStatus\n");
  }

  /**
   * Write the new alleles section to the supplied output
   * @param output the stream to write to
   * @param alleleTranslation the allele translation to read from
   * @throws IOException
   */
  private void writeAlleles(BufferedWriter output, AlleleTranslation alleleTranslation) throws IOException {
    String namedAlleles = alleleTranslation.getNamedAlleles().stream()
        .map(a -> {
          String alleles = Arrays.stream(a.getAlleles())
              .map(b -> b == null ? "" : b).collect(Collectors.joining("\t"));
          return "Allele\t" + a.getId() + "\t" + a.getName() + "\t" + a.getFunction() + "\t" + alleles;
        })
        .collect(Collectors.joining("\n"));
    output.write(namedAlleles);
    output.write("\n");
  }

  /**
   * Write the new notes section to the supplied output
   * @param output the stream to write to
   * @param alleleTranslation the allele translation to read from
   * @throws IOException
   */
  private void writeNotes(BufferedWriter output, AlleleTranslation alleleTranslation) throws IOException {
    if (alleleTranslation.getNotes() == null) {
      return;
    }

    String notes = alleleTranslation.getNotes().stream()
        .map(n -> "Note\t" + n + "\n")
        .collect(Collectors.joining());
    output.write(notes);
    output.write("\n");
  }
}
