package org.pharmgkb.pharmcat.definition;

import java.io.InputStream;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.net.HttpURLConnection;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.definition.model.VariantType;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.Iupac;
import org.pharmgkb.pharmcat.util.DataManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.w3c.dom.Document;


/**
 * This extracts the positions from the definition files and reformats them as a vcf. On the command line it takes one
 * argument:
 *
 * <ul>
 *   <li>-o output_vcf_path = A path to write the VCF file to</li>
 * </ul>
 *
 * @author Lester Carter
 */
public class ExtractPositions {

  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final Path sf_definitionDir = DataManager.DEFAULT_DEFINITION_DIR;
  private static final String sf_fileHeader = "##fileformat=VCFv4.1\n" +
      "##fileDate=%s\n" +
      "##source=PharmCAT allele definitions\n" +
      "##reference=hg38\n" +
      "##contig=<ID=chr1,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr2,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr3,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr4,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr5,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr6,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr7,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr8,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr9,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr10,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr11,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr12,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr13,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr14,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr15,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr16,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr17,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr18,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr19,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr20,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr21,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chr22,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chrX,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##contig=<ID=chrY,assembly=hg38,species=\"Homo sapiens\">\n" +
      "##INFO=<ID=PX,Number=.,Type=String,Description=\"PGX\">\n" +
      "##INFO=<ID=POI,Number=0,Type=Flag,Description=\"Position of Interest but not part of an allele definition\">\n" +
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
      "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n";

  private final Path m_outputVcf;
  private final DocumentBuilderFactory f_dbf;


  // Default constructor
  public ExtractPositions(Path outputVcf) {
    m_outputVcf = outputVcf;
    f_dbf = DocumentBuilderFactory.newInstance();
    f_dbf.setValidating(false);
    try {
      f_dbf.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false);
    } catch (ParserConfigurationException e) {
      throw new RuntimeException("Error setting up XML parser", e);
    }
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

      ExtractPositions extractPositions = new ExtractPositions(outputVcf);
      extractPositions.run();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  //Extract the positions and print them to file
  private void run() {
    try {
      DefinitionReader definitionReader = new DefinitionReader();
      definitionReader.read(sf_definitionDir);
      if (definitionReader.getGenes().size() == 0) {
        sf_logger.error("Did not find any allele definitions at {}", sf_definitionDir);
        System.exit(1);
      }
      sf_logger.info("Writing to {}", m_outputVcf);
      StringBuilder vcfString = sortVcf(getPositions(definitionReader));
      try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(m_outputVcf))) {
        writer.print(vcfString);
        writer.flush();
      }
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }

  }


  // Build up vcf string
  public StringBuilder getPositions(@Nonnull DefinitionReader definitionReader) {
    StringBuilder builder = new StringBuilder();
    builder.append(String.format(sf_fileHeader, DateTimeFormatter.ISO_OFFSET_DATE_TIME.format(ZonedDateTime.now())));
    for (String gene : definitionReader.getGenes()) {  // for each definition file
      int positionCount = 0;
      // convert bXX format to hgXX
      String build = "hg" + definitionReader.getDefinitionFile(gene).getGenomeBuild().substring(1);
      for (VariantLocus variantLocus : definitionReader.getPositions(gene)) {
        positionCount++;
        String[] vcfFields = getVcfLineFromDefinition(definitionReader, gene, variantLocus, build);
        String vcfLine = String.join("\t", vcfFields);
        builder.append(vcfLine).append("\n");
      }
      DefinitionExemption exemption = definitionReader.getExemption(gene);
      if (exemption != null) {
        for (VariantLocus variant : exemption.getExtraPositions()) {
          positionCount++;
          String[] vcfFields = getVcfLineFromDefinition(definitionReader, gene, variant, build);
          String vcfLine = String.join("\t", vcfFields);
          builder.append(vcfLine).append("\n");
        }
      }
      sf_logger.info("queried {} positions for {}", positionCount, gene);
    }
    return  builder;
  }


  /**
   * Helper method to convert repeats into standard format
   */
  private static String expandRepeats(@Nonnull String allele) {
    String finalAllele = allele;
    if (allele.contains("(")) {
      int bracketStart=0;
      int bracketEnd=0;
      StringBuilder repeat = new StringBuilder();
      for (int i = 0; i < allele.toCharArray().length; i++) {
        String character = String.valueOf(allele.toCharArray()[i]);
        if (character.equals("(")) {
          bracketStart = i;
        }
        if (character.equals(")")) {
          bracketEnd = i;
        }
        if (character.matches("[0-9]+")) {
          repeat.append(character); // just in case we have a repeat greater than nine

        }
      }
      String repeatSeq = allele.substring(bracketStart+1, bracketEnd);
      StringBuilder expandedAlelle = new StringBuilder();
      for (int repeats = 0; repeats< Integer.parseInt(repeat.toString()); repeats++) {
        expandedAlelle.append(repeatSeq);
      }
      finalAllele = allele.substring(0, bracketStart) + expandedAlelle + allele.substring(bracketEnd+1+repeat.length());
    }
    return finalAllele;
  }


  /**
   * For deletions we need to be able to get the nucleotide at the previous position.
   * This is not stored in the definition file, so we need to use an external source.
   */
  String getDAS(String chr, String position, String genomeBuild) {
    String nucleotide;
    try {
      String uri =
          "http://genome.ucsc.edu/cgi-bin/das/"+ genomeBuild +  "/dna?segment="+ chr +":" +position +"," + position;
      sf_logger.debug("Getting position: {}", uri);
      URL url = new URL(uri);
      HttpURLConnection connection =
          (HttpURLConnection) url.openConnection();
      connection.setRequestMethod("GET");
      connection.setRequestProperty("Accept", "application/xml");
      InputStream xml = connection.getInputStream();
      DocumentBuilder db = f_dbf.newDocumentBuilder();
      Document doc = db.parse(xml);
      doc.getDocumentElement().normalize();
      nucleotide = doc.getElementsByTagName("DNA").item(0).getTextContent();
    } catch (Exception e) {
      throw new RuntimeException("Error when requesting DAS data", e);
    }
    nucleotide = nucleotide.toUpperCase().trim();
    return nucleotide;
  }


  /**
   * Simple sorting of Vcf file on chr and position
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





  /**
   * Helper method to convert repeats into standard format
   */
  private String[] getVcfLineFromDefinition(@Nonnull DefinitionReader definitionReader, @Nonnull String gene,
      @Nonnull VariantLocus variantLocus, String genomeBuild) {

    DefinitionFile definitionFile = definitionReader.getDefinitionFile(gene);
    NamedAllele namedAlleleFirst = definitionFile.getNamedAlleles().get(0); //Get star one for Ref column of vcf. Use DAS in future
    String allele = "";

    if (namedAlleleFirst.getAllele(variantLocus) != null) {  // for the exceptions/exemptions
      allele = namedAlleleFirst.getAllele(variantLocus);
    }
    int position = variantLocus.getVcfPosition();
    String rsid = variantLocus.getRsid();
    if (rsid == null) {
      rsid = ".";
    }
    String chr = definitionReader.getDefinitionFile(gene).getChromosome();
    List<String> alts = new ArrayList<>(); //get alts from the namedAlleles

    String dasRef = getDAS(chr, String.valueOf(position), genomeBuild);
    // allele may be blank if the position is an "extra" position not used in an allele definition
    if (StringUtils.isBlank(allele)) {
      allele = dasRef;
    }
    
    if (!dasRef.equals(allele) && variantLocus.getType()== VariantType.SNP && allele.length()==1) {
      sf_logger.warn("Star one/ref mismatch at position: {} - using {} as ref as opposed to {}", position, dasRef, allele);
      alts.add(allele);
      allele = dasRef;
    }

    List<String> starAlleles = new ArrayList<>();
    String finalAlleleRef = allele;
    definitionFile.getNamedAlleles().stream()
        .filter(namedAllele -> namedAllele.getAllele(variantLocus) != null)
        .forEachOrdered(namedAllele -> {
          if (!alts.contains(expandRepeats(namedAllele.getAllele(variantLocus)))
              && !finalAlleleRef.equals(expandRepeats(namedAllele.getAllele(variantLocus)))) {
            alts.add(expandRepeats(namedAllele.getAllele(variantLocus)));
          }

      //calculate size of haplotype
      List<String> l = Arrays.stream(namedAllele.getAlleles())
          .filter(Objects::nonNull)
          .collect(Collectors.toList());
      starAlleles.add(gene + ":" + namedAllele.getName().replace(" ","")+"[" + l.size() +"]" +"is" + namedAllele.getAllele(variantLocus));
    });
    if (alts.size() == 0 ) {
      alts.add(".");
    }
    String finalAllele = expandRepeats(allele);
    // Simplest possible del/ins parsing, presuming the only a single ins or del string, or the correct one being first
    if (alts.size()>0) {
      if (alts.stream().anyMatch(a -> a.contains("ins"))) {
        String nucleotide = getDAS(chr, Integer.toString(position), genomeBuild);
        for (int i = 0; i < alts.size(); i++) {
          alts.set(i, alts.get(i).replace("ins", nucleotide));
        }
        finalAllele = nucleotide;
      }
      if (alts.get(0).contains("del")) {
        String nucleotide = getDAS(chr, Integer.toString(position), genomeBuild);
        finalAllele = nucleotide+ alts.get(0).replace("del", "");
        alts.set(0, nucleotide);
      }
    }
    SortedSet<String> expandedAlts = expandIupacAlts(alts);
    alts.remove(finalAllele);

    String alt = String.join(",", expandedAlts);
    String starAllele;
    if (starAlleles.size() > 0) {
      starAllele = "PX=" + String.join(",", starAlleles);
    } else {
      starAllele = "POI";
    }

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

  /**
   * For the list of given alt alleles, replace Iupac ambiguity codes with all the bases they represent
   * @param alts a Collection of alt alleles
   * @return a SortedSet of alt alleles with no ambiguity codes (may include the ref allele)
   */
  private static SortedSet<String> expandIupacAlts(Collection<String> alts) {
    SortedSet<String> bases = new TreeSet<>();

    for (String alt : alts) {
      if (alt.length() == 1 && !alt.equals(".")) {
        try {
          Iupac iupac = Iupac.lookup(alt);
          bases.addAll(iupac.getBases());
        } catch (IllegalArgumentException ex) {
          throw new RuntimeException("Bad IUPAC code " + alt, ex);
        }
      } else {
        bases.add(alt);
      }
    }

    return bases;
  }
}
