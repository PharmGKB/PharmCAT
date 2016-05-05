package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;


/**
 * Serializer/Deserializer for final allele definition file.
 *
 * @author Mark Woon
 */
public class GeneratedDefinitionSerializer {
  private static final Gson sf_gson = new GsonBuilder().serializeNulls().setPrettyPrinting().create();
  private SimpleDateFormat m_dateFormat = new SimpleDateFormat("MM/dd/yy");


  public void serializeToJson(@Nonnull DefinitionFile definitionFile, @Nonnull Path jsonFile) throws IOException {

    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"), "Invalid format: file name does not end with .json");

    try (BufferedWriter writer = Files.newBufferedWriter(jsonFile, StandardCharsets.UTF_8)) {
      sf_gson.toJson(definitionFile, writer);
    }
  }


  public DefinitionFile deserializeFromJson(@Nonnull Path jsonFile) throws IOException {

    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"), "Invalid format: file name does not end with .json");

    try (BufferedReader reader = Files.newBufferedReader(jsonFile, StandardCharsets.UTF_8)) {
      return sf_gson.fromJson(reader, DefinitionFile.class);
    }
  }


  public void serializeToTsv(@Nonnull DefinitionFile definitionFile, @Nonnull Path tsvFile) throws IOException {

    Preconditions.checkNotNull(tsvFile);
    Preconditions.checkArgument(tsvFile.toString().endsWith(".tsv"), "Invalid format: file name does not end with .tsv");

    try (PrintWriter output = new PrintWriter(Files.newBufferedWriter(tsvFile, StandardCharsets.UTF_8))) {

      writeHeader(output, definitionFile);
      writeAlleles(output, definitionFile);
      writeNotes(output, definitionFile);

    } catch (IOException e) {
      throw new RuntimeException("Couldn't write output", e);
    }
  }

  /**
   * Write the new header section to the supplied output.
   * @param output the stream to write to
   * @param definitionFile the allele translation to read from
   */
  private void writeHeader(PrintWriter output, DefinitionFile definitionFile) throws IOException {
    output.println("FormatVersion\t"    + definitionFile.getFormatVersion());
    output.println("GeneName\t"         + definitionFile.getGeneSymbol());
    output.println("GeneRefSeq\t"       + definitionFile.getRefSeqGene());
    output.println("GeneOrientation\t"  + definitionFile.getOrientation());
    output.println("ContentDate\t"      + m_dateFormat.format(definitionFile.getModificationDate()));
    output.println("ContentVersion\t"   + definitionFile.getContentVersion());
    output.println("GenomeBuild\t"      + definitionFile.getGenomeBuild());
    output.println("ChrName\t"          + definitionFile.getChromosome());
    output.println("ChrRefSeq\t"        + definitionFile.getRefSeqChromosome());
    output.println("ProteinRefSeq\t"    + definitionFile.getRefSeqProtein());
    output.println("NumVariants\t"      + definitionFile.getVariants().length);

    output.print("ResourceNote\t\t\t\t");
    output.println(Arrays.stream(definitionFile.getVariants()).map(VariantLocus::getResourceNote).collect(Collectors.joining("\t")));

    output.print("ProteinNote\t\t\t\t");
    output.println(Arrays.stream(definitionFile.getVariants()).map(VariantLocus::getProteinNote).collect(Collectors.joining("\t")));

    output.print("ChrPosition\t\t\t\t");
    output.println(Arrays.stream(definitionFile.getVariants()).map(VariantLocus::getChromosomeHgvsName).collect(Collectors.joining("\t")));

    output.print("GenePosition\t\t\t\t");
    output.println(Arrays.stream(definitionFile.getVariants()).map(VariantLocus::getGeneHgvsName).collect(Collectors.joining("\t")));

    output.print("rsID\t\t\t\t");
    output.println(Arrays.stream(definitionFile.getVariants()).map(VariantLocus::getRsid).collect(Collectors.joining("\t")));

    output.println("Header\tID\tName\tFunctional Status");
  }


  private static final Joiner sf_alleleJoiner = Joiner.on("\t").useForNull("");

  /**
   * Write the new alleles section to the supplied output.
   *
   * @param output the stream to write to
   * @param definitionFile the allele translation to read from
   */
  private void writeAlleles(PrintWriter output, DefinitionFile definitionFile) throws IOException {

    for (NamedAllele allele : definitionFile.getNamedAlleles()) {
      output.print("Allele\t");
      output.print(allele.getId());
      output.print("\t");
      output.print(allele.getName());
      output.print("\t");
      output.print(allele.getFunction());
      output.print("\t");
      output.println(sf_alleleJoiner.join(allele.getAlleles()));
    }
    output.println();
  }

  /**
   * Write the new notes section to the supplied output.
   *
   * @param output the stream to write to
   * @param definitionFile the allele translation to read from
   */
  private void writeNotes(PrintWriter output, DefinitionFile definitionFile) throws IOException {

    if (definitionFile.getNotes() == null) {
      return;
    }
    for (String note : definitionFile.getNotes()) {
      output.print("Note\t");
      output.println(note);
    }
    output.println();
  }
}
