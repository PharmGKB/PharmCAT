package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import javax.annotation.Nonnull;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.util.DataSerializer;


/**
 * Created by lester on 1/2/17.
 *
 * Feasibility study of liftover on the definition files.
 *
 * Still very much experimental work in progress.  For use with the grc37 test data sets.  Talk to you Ryan for details.
 *
 * User runs to get a grc38 bed, and uses this as in input for remap using
 * ncbi remap: https://www.ncbi.nlm.nih.gov/genome/tools/remap
 *
 * Remap using homo sapien GRCh38.p6 to GRCh37.p13.  Use the bed file input and output, and use defaults for the rest.
 *
 * The user downloads the resulting grc37 bed and it is fed back into the program to make the  grc37 definitions.
 *
 */


public class ConvertDefinition {
  private static Path m_definitionDir;
  private static DefinitionReader m_definitionReader;
  private static Path m_outputDefinitionDir;
  private static String m_genomeBuildFrom;
  private static String m_genomeBuildTo;

  // Default constructor
  public ConvertDefinition (Path definitionDir, Path outputDefinitionDir, String from, String to) {
    m_definitionDir = definitionDir;
    m_outputDefinitionDir =outputDefinitionDir;
    m_genomeBuildFrom = from;
    m_genomeBuildTo = to;
  }


  // Main method for command line use
  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("d", "definition-dir", "directory of allele definition grc38 json files", true, "d")
          .addOption("o", "output-dif", "output directory of allele definition grc37 json files", true, "o")
          .addOption("f", "from", "from", true, "f")
          .addOption("t", "to", "to", true, "t");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }
      Path definitionDir = cliHelper.getValidDirectory("d", false);
      Path outputDefinitionDir= cliHelper.getValidDirectory("o", false);
      ConvertDefinition convertDefintion = new ConvertDefinition(definitionDir, outputDefinitionDir,
          cliHelper.getValue("f"),  cliHelper.getValue("t"));
      convertDefintion.run();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  private void run() {
    System.out.println("Running convert definitions");
    try {
      m_definitionReader = new DefinitionReader();
      System.out.println(m_definitionDir);
      m_definitionReader.read(m_definitionDir);
      if (m_definitionReader.getGenes().size() == 0) {
        System.out.println("Did not find any allele definitions at " + m_definitionDir);
      }
      StringBuilder intervalString = getPositions(m_definitionReader);

      // todo: make Win etc safe
      System.out.println(m_outputDefinitionDir.toString()+"/temp.bed");
      try (PrintWriter writer = new PrintWriter(m_outputDefinitionDir.toString()+"/temp.bed",
          "UTF-8")) {  // close PrintWriter with try
        writer.print(intervalString);
        writer.flush();
      }

      Scanner s = new Scanner(System.in);
      System.out.println("A .bed file has been created in the output directory.  Please run remap and " +
          "save as remap.bed. Press return when done");
      s.nextLine();
      Path remapBed = Paths.get(m_outputDefinitionDir.toString()+"/remap.bed");

      while (! Files.exists(remapBed)) {
        System.out.println("Not valid, try again");
        s.nextLine();
      }
      updateDefinitions(remapBed);
      writeNewDefinitions();
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  private void writeNewDefinitions() {
    for (String gene : m_definitionReader.getGenes()) {
      DefinitionFile definitionFile = m_definitionReader.getDefinitionFile(gene);
      Path jsonFile = Paths.get(m_outputDefinitionDir+ "/" + gene +"_translation.json");
      DataSerializer serializer = new DataSerializer();
      try {
        serializer.serializeToJson(definitionFile, jsonFile);
      } catch (IOException e) {
        e.printStackTrace();
      }
      System.out.println("Wrote " + jsonFile);
      System.out.println();
    }
  }

  private void updateDefinitions(Path remapBed) {
    try {
      BufferedReader in = new BufferedReader(new FileReader(remapBed.toString()));
      String line;
      while ((line = in.readLine()) != null) {
        ArrayList<String> lineDetails = new ArrayList<>(Arrays.asList(line.split("\t")));
        String oldPosition = lineDetails.get(3).split("_")[1];
        String newPosition = lineDetails.get(1);
        String gene = lineDetails.get(3).split("_")[2];


        for (VariantLocus variantLocus : m_definitionReader.getPositions(gene)) {
          if (Integer.parseInt(oldPosition) == variantLocus.getPosition()) {
            System.out.println("old->new ("+ gene + "): " +oldPosition + "\t" + newPosition);
            variantLocus.setPosition(Integer.parseInt(newPosition));
            variantLocus.setResourceNote(variantLocus.getResourceNote()+ " Updated position "
                + "old->new ("+ gene + "): " +oldPosition + " to " + newPosition);
          }
        }
      }
      in.close();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }


  // Build up string
  public  StringBuilder getPositions(@Nonnull DefinitionReader definitionReader) throws IOException {
    StringBuilder builder = new StringBuilder();
    for (String gene : definitionReader.getGenes()) {  // for each definition file
      for (VariantLocus variantLocus :definitionReader.getPositions(gene)) {
        int start = variantLocus.getPosition();
        int end = variantLocus.getPosition()+1;
        String summary_info = definitionReader.getDefinitionFile(gene).getChromosome() +"_" + start+ "_" + gene;
        String interval = definitionReader.getDefinitionFile(gene).getChromosome() + "\t" +
            start + "\t" + end + "\t" + summary_info;
        System.out.println("Gene:" + gene + " - " + interval);
        builder.append(interval).append("\n");
      }
    }
    return  builder;
  }


}

