package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import java.util.TreeMap;
import javax.annotation.Nonnull;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.definition.model.VariantType;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.util.DataManager;
import org.w3c.dom.Document;


/**
 * This extracts the positions from the definition files and reformats them as a vcf.
 * As a command line it takes three arguments:
 * -d  CPIC definitions directory
 * -o output vcf
 * -g genome build to use for DAS
 *
 * @author lester
 *
 */

public class ExtractPositions {
  private static Path ps_definitionDir;
  private static Path ps_outputVcf;
  private static final String sf_fileHeader = "##fileformat=VCFv4.1\n" +
      "##fileDate=2015-08-04\n" +
      "##source=Electronic, version: hg38_2.0.1\n" +
      "##reference=hg38\n" +
      "##INFO=<ID=PX,Number=.,Type=String,Description=\"PGX\">\n" +
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
      "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n";


  // Default constructor
  public ExtractPositions(Path definitionDir, Path outputVcf) {
    ps_definitionDir = definitionDir;
    ps_outputVcf = outputVcf;
  }


  // Main method for command line use
  public static void main(String[] args) {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("o", "output-file", "output vcf file", true, "o");
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }
      Path outputVcf= cliHelper.getValidFile("o", false);

      ExtractPositions extractPositions = new ExtractPositions(DataManager.DEFAULT_DEFINITION_DIR, outputVcf);
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
      StringBuilder vcfString = sortVcf(getPositions(definitionReader));
      try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(ps_outputVcf))) {
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
      // convert bXX format to hgXX
      String build = "hg" + definitionReader.getDefinitionFile(gene).getGenomeBuild().substring(1);
      for (VariantLocus variantLocus : definitionReader.getPositions(gene)) {
        String[] vcfFields = getVcfLineFromDefinition(definitionReader, gene, variantLocus, build);
        String vcfLine = String.join("\t", (CharSequence[])vcfFields);
        builder.append(vcfLine).append("\n");
      }
      DefinitionExemption exemption = definitionReader.getExemption(gene);
      for (String rsid : exemption.getExtraPositions()) {

      }
    }
    return  builder;
  }


  /*
  * Helper method to convert repeats into standard format
  */
  public static String expandAllele(@Nonnull String allele) {
    String finalAllele = allele;
    if (allele.contains("(")) {
      int bracketStart=0;
      int bracketEnd=0;
      String repeat = "";
      for (int i = 0; i < allele.toCharArray().length; i++) {
        String character = String.valueOf(allele.toCharArray()[i]);
        if (character.equals("(")) {
          bracketStart = i;
        }
        if (character.equals(")")) {
          bracketEnd = i;
        }
        if (character.matches("[0-9]+")) {
          repeat = repeat + character; // just in case we have a repeat greater than nine

        }
      }
      String repeatSeq = allele.substring(bracketStart+1, bracketEnd);
      String expandedAlelle= "";
      for (int repeats =0; repeats< Integer.parseInt(repeat); repeats++) {
        expandedAlelle = expandedAlelle + repeatSeq;
      }
      finalAllele = allele.substring(0, bracketStart) + expandedAlelle + allele.substring(bracketEnd+1+repeat.length(), allele.length());
    }
    return finalAllele;
  }


  /*
  * For deletions we need to be able to get the nucleotide at the previous position.
  * This is not stored in the definition file, so we need to use an external source.
  *
   */
  public static String getDAS(String chr, String position, String genomeBuild) {
    String nucleotide= "";
    try {
      String uri =
          "http://genome.ucsc.edu/cgi-bin/das/"+ genomeBuild +  "/dna?segment="+ chr +":" +position +"," + position;
      System.out.println("Getting position: " + uri);
      URL url = new URL(uri);
      HttpURLConnection connection =
          (HttpURLConnection) url.openConnection();
      connection.setRequestMethod("GET");
      connection.setRequestProperty("Accept", "application/xml");
      InputStream xml = connection.getInputStream();
      DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
      DocumentBuilder db = dbf.newDocumentBuilder();
      Document doc = db.parse(xml);
      doc.getDocumentElement().normalize();
      nucleotide = doc.getElementsByTagName("DNA").item(0).getTextContent();
    }  catch (Exception e) {
      e.printStackTrace();
    }
    nucleotide = nucleotide.toUpperCase().trim();
    return nucleotide;
  }


  /*
  * Simple sorting of Vcf file on chr and position
  *
   */
  private StringBuilder sortVcf (StringBuilder vcfString) {
    StringBuilder sortedVcfString = new StringBuilder();
    Map<Integer, Map<Integer, String>> treeMapChr = new TreeMap<>();
    for (String vcfLine: vcfString.toString().split("\n")) {
      if (!vcfLine.startsWith("#")){
        String[] vcfFields = vcfLine.split("\t");
        Map<Integer, String> treeMapPosition = new TreeMap<>();
        int chr = Integer.parseInt(vcfFields[0].replace("chr", ""));
        int position = Integer.parseInt(vcfFields[1]);
        if (treeMapChr.containsKey(chr)) {
          Map<Integer, String> treeMapPositionCurrent = treeMapChr.get(chr);
          treeMapPositionCurrent.put(position, vcfLine);
          treeMapChr.put(chr, treeMapPositionCurrent);
        }
        else {
          treeMapPosition.put(position, vcfLine);
          treeMapChr.put(chr, treeMapPosition);
        }
      }
      else {
        sortedVcfString.append(vcfLine).append("\n");
      }
    }

    for (Map<Integer, String> chrMap: treeMapChr.values()) {
      for (Object posMap: chrMap.values()) {
        sortedVcfString.append(posMap.toString()).append("\n");
      }
    }
    return sortedVcfString;
  }





  /*
  * Helper method to convert repeats into standard format
  */
  public static String[] getVcfLineFromDefinition(@Nonnull DefinitionReader definitionReader, @Nonnull String gene, @Nonnull VariantLocus variantLocus,
      String genomeBuild) {
    DefinitionFile definitionFile = definitionReader.getDefinitionFile(gene);
    NamedAllele namedAlleleFirst = definitionFile.getNamedAlleles().get(0); //Get star one for Ref column of vcf. Use DAS in future
    String allele = namedAlleleFirst.getAllele(variantLocus);
    int position = variantLocus.getVcfPosition();
    String rsid = variantLocus.getRsid();
    if (rsid == null) {
      rsid = ".";
    }
    String chr = definitionReader.getDefinitionFile(gene).getChromosome();
    ArrayList<String> alts = new ArrayList<>(); //get alts from the namedAlleles
    String dasRef = getDAS(chr, String.valueOf(position), genomeBuild);
    if (!dasRef.equals(allele) && variantLocus.getType()== VariantType.SNP && allele.length()==1) {
      System.out.println("Star one/ref mismatch at position: " + position + " - using " + dasRef + " as ref as opposed to " + allele);
      alts.add(allele);
      allele = dasRef;
    }
    ArrayList<String> starAlleles = new ArrayList<>();
    String finalAlleleRef = allele;
    definitionFile.getNamedAlleles().stream().filter(namedAllele -> namedAllele.getAllele(variantLocus) != null
    ).forEachOrdered(namedAllele -> {
      if (!alts.contains(expandAllele(namedAllele.getAllele(variantLocus)))
          && !finalAlleleRef.equals(expandAllele(namedAllele.getAllele(variantLocus)))) {
        alts.add(expandAllele(namedAllele.getAllele(variantLocus)));
      }

      //calculate size of haplotype
      ArrayList<String> l = new ArrayList<>(Arrays.asList(namedAllele.getAlleles()));
      l.removeAll(Collections.singleton(null));
      starAlleles.add(gene + ":" + namedAllele.getName().replace(" ","")+"[" + l.size() +"]" +"is" + namedAllele.getAllele(variantLocus));
    });
    if (alts.size() == 0 ) {
      alts.add(".");
    }
    String finalAllele = expandAllele(allele);
    alts.remove(finalAllele); // get rid of any lingering duplicates
    if (alts.contains("Y")) {
      alts.remove("Y");
      if (!alts.contains("C") && !allele.equals("C")) {
        alts.add("C");
      }
      if (!alts.contains("T")  && !allele.equals("T")) {
        alts.add("T");
      }
    }
    if (alts.contains("R")) {
      alts.remove("R");
      if (!alts.contains("G")  && !allele.equals("G")) {
        alts.add("G");
      }
      if (!alts.contains("A")  && !allele.equals("A")) {
        alts.add("A");
      }
    }
    // Simplest possible del/ins parsing, presuming the only a single ins or del string, or the correct one being first
    if (alts.size()>0) {
      if (alts.get(0).contains("ins")) {
        String nucleotide = getDAS(chr, Integer.toString(position), genomeBuild);
        alts.set(0, nucleotide + alts.get(0).replace("ins", ""));
        finalAllele = nucleotide;
      }
      if (alts.get(0).contains("del")) {
        String nucleotide = getDAS(chr, Integer.toString(position), genomeBuild);
        finalAllele = nucleotide+ alts.get(0).replace("del", "");
        alts.set(0, nucleotide);
      }
    }
    String alt = String.join(",", alts);
    String starAllele = "PX=" + String.join(",", starAlleles)+";";

    return new String[]{
        chr,
        Integer.toString(position),
        rsid,
        finalAllele,
        alt,
        ".",
        "PASS",
        starAllele,
        "GT",
        "0/0" };
  }
}
