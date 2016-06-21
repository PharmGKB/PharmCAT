package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.ArrayList;
import javax.annotation.Nonnull;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;


/**
 * This extracts the positions from the definition files and reformats them as a vcf.
 * As a command line it takes two arguments:
 * -d  CPIC definitions directory
 * -o output vcf
 *
 * Created by lester on 6/20/16.
 */

public class ExtractPositions {
  private static Path ps_definitionDir;
  private static Path ps_outputVcf;
  private static final String sf_fileHeader = "##fileformat=VCFv4.1\n" +
      "##fileDate=2015-08-04\n" +
      "##source=Electronic, version: hg38_2.0.1\n" +
      "##reference=hg38\n" +
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
      "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
      "#CHROM  POS  ID REF  ALT  QUAL FILTER INFO FORMAT PharmCAT\n";


  // Default constructor
  public ExtractPositions(Path definitionDir, Path outputVcf) {
    ps_definitionDir = definitionDir;
    ps_outputVcf = outputVcf;
  }


  // Main method for command line use
  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("d", "definition-dir", "directory of allele definition files", true, "d")
          .addOption("o", "output-file", "output vcf file", true, "o");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }
      Path definitionDir = cliHelper.getValidDirectory("d", false);
      Path outputVcf= cliHelper.getValidFile("o", false);
      ExtractPositions extractPositions = new ExtractPositions(definitionDir, outputVcf);
      extractPositions.run();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  //Extract the positions and print them to file
  private void run() {
    try {
      DefinitionReader definitionReader = new DefinitionReader();
      definitionReader.read(ps_definitionDir);
      if (definitionReader.getGenes().size() == 0) {
        System.out.println("Did not find any allele definitions at " + ps_definitionDir);
        System.exit(1);
      }
      StringBuilder vcfString = getPositions(definitionReader);
      try (PrintWriter writer = new PrintWriter(String.valueOf(ps_outputVcf), "UTF-8")) {  // close PrintWriter with try
        writer.print(vcfString);
        writer.flush();
      }
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }

  }


  // Build up vcf string
  public  StringBuilder getPositions(@Nonnull DefinitionReader definitionReader) throws IOException {
    StringBuilder builder = new StringBuilder();
    builder.append(sf_fileHeader);
    for (String gene : definitionReader.getGenes()) {  // for each definition file
      DefinitionFile definitionFile = definitionReader.getDefinitionFile(gene);
      VariantLocus[] positions = definitionReader.getPositions(gene);
      NamedAllele namedAlleleFirst = definitionFile.getNamedAlleles().get(0); //Get star one for Ref column of vcf
      int countPosition = 0;
      for (String allele :namedAlleleFirst.getAlleles() ) {
        if (allele != null) {
          int position = positions[countPosition].getVcfPosition();
          String rsid;
          rsid = positions[countPosition].getRsid();
          if (rsid == null) {
            rsid = ".";
          }
          ArrayList <String> alts = new ArrayList<>(); //get alts from the namedAlleles
          ArrayList <String> starAlleles = new ArrayList<>();
          for (NamedAllele namedAllele: definitionFile.getNamedAlleles()) {
            if (namedAllele.getAlleles()[countPosition] != null && !namedAllele.getAlleles()[countPosition].equals(allele)
                && !alts.contains(namedAllele.getAlleles()[countPosition])) {
              alts.add(namedAllele.getAlleles()[countPosition]);
              starAlleles.add(gene + ":" +namedAllele.getName());
            }
          }
          String alt = String.join(",", alts);
          String starAllele = String.join(",", starAlleles);

          String[] vcfFields = {
              definitionFile.getChromosome(),
              Integer.toString(position),
              rsid,
              allele,
              alt,
              ".",
              "PASS",
              starAllele,
              "GT",
              "0/0\n",
          };
          String vcfLine = String.join("\t", (CharSequence[])vcfFields);
          builder.append(vcfLine);
        }
        countPosition = countPosition+1;
      }
    }
    return  builder;
  }
}

