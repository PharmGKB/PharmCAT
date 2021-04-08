package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
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
 * @author Ryan Whaley
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
      "##INFO=<ID=PX,Number=.,Type=String,Description=\"Gene(:Allele Name(n of defining positions)defining allele)\">\n" +
      "##INFO=<ID=POI,Number=0,Type=Flag,Description=\"Position of Interest but not part of an allele definition\">\n" +
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
      "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n";

  private final Path m_outputVcf;
  private final DocumentBuilderFactory f_dbf;
  private final DefinitionReader f_definitionReader;
  private final Map<String,String> f_refCache = new HashMap<>();


  /**
   * Default constructor
   * @param outputVcf file path to write output VCF to
   */
  public ExtractPositions(Path outputVcf) {
    m_outputVcf = outputVcf;
    f_dbf = DocumentBuilderFactory.newInstance();
    f_dbf.setValidating(false);
    try {
      f_dbf.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false);

      // load the allele definition files
      f_definitionReader = new DefinitionReader();
      f_definitionReader.read(sf_definitionDir);
      if (f_definitionReader.getGenes().size() == 0) {
        throw new RuntimeException("Did not find any allele definitions at " + sf_definitionDir);
      }

      // load the DAS lookup cache
      try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("positions_reference.tsv")))) {
        String line;
        List<String> errors = new ArrayList<>();
        while ((line = reader.readLine()) != null) {
          String[] fields = line.split("\t");
          String previousValue = f_refCache.put(fields[0], fields[1]);
          if (previousValue != null) {
            errors.add("Duplicate key found: " + fields[0]);
          }
        }
        if (errors.size() > 0) {
          throw new RuntimeException("Found problems in the reference cache\n" + String.join("\n", errors));
        }
      }
    } catch (ParserConfigurationException e) {
      throw new RuntimeException("Error setting up XML parser", e);
    } catch (IOException e) {
      throw new RuntimeException("Error reading gene definitions", e);
    }
  }


  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("o", "output-file", "output vcf file", true, "o");

    cliHelper.execute(args, cli -> {
      try {
        ExtractPositions extractPositions = new ExtractPositions(cli.getValidFile("o", false));
        extractPositions.run();
        return 0;
      } catch (Exception ex) {
        ex.printStackTrace();
        return 1;
      }
    });
  }


  /**
   * Extract the positions and print them to file.
   */
  private void run() {
    try {
      StringBuilder vcfString = sortVcf(getPositions());
      try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(m_outputVcf))) {
        sf_logger.info("Writing to {}", m_outputVcf);
        writer.print(String.format(sf_fileHeader, DateTimeFormatter.ISO_OFFSET_DATE_TIME.format(ZonedDateTime.now())));
        writer.print(vcfString);
        writer.flush();
      }
      sf_logger.info("Done");
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  /**
   * Build up VCF positions.
   */
  public StringBuilder getPositions() {
    StringBuilder builder = new StringBuilder();
    for (String gene : f_definitionReader.getGenes()) {  // for each definition file
      DefinitionFile definitionFile = f_definitionReader.getDefinitionFile(gene);
      int positionCount = 0;
      // convert bXX format to hgXX
      String build = "hg" + definitionFile.getGenomeBuild().substring(1);
      for (VariantLocus variantLocus : definitionFile.getVariants()) {
        positionCount++;
        String[] vcfFields = getVcfLineFromDefinition(definitionFile, variantLocus, build);
        String vcfLine = String.join("\t", vcfFields);
        builder.append(vcfLine).append("\n");
      }
      DefinitionExemption exemption = f_definitionReader.getExemption(gene);
      if (exemption != null) {
        for (VariantLocus variant : exemption.getExtraPositions()) {
          positionCount++;
          String[] vcfFields = getVcfLineFromDefinition(definitionFile, variant, build);
          String vcfLine = String.join("\t", vcfFields);
          builder.append(vcfLine).append("\n");
        }
      }
      sf_logger.info("queried {} positions for {}", positionCount, gene);
    }
    return builder;
  }


  /**
   * Helper method to convert repeats into standard format
   */
  private static String expandRepeats(String allele) {
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

    // check the cache before doing an expensive DAS lookup
    String cacheKey = genomeBuild + ":" + chr + ":" + position;
    if (f_refCache.containsKey(cacheKey)) {
      return f_refCache.get(cacheKey);
    }

    // do the call out to UCSC to get the ref sequence
    try {
      String uri =
          "http://genome.ucsc.edu/cgi-bin/das/" + genomeBuild + "/dna?segment=" + chr + ":" + position + "," + position;
      sf_logger.debug("Getting position: {}", uri);
      URL url = new URL(uri);
      HttpURLConnection connection = (HttpURLConnection) url.openConnection();
      connection.setRequestMethod("GET");
      connection.setRequestProperty("Accept", "application/xml");
      InputStream xml = connection.getInputStream();
      DocumentBuilder db = f_dbf.newDocumentBuilder();
      Document doc = db.parse(xml);
      doc.getDocumentElement().normalize();
      String nucleotide = doc.getElementsByTagName("DNA").item(0).getTextContent();
      nucleotide = nucleotide.toUpperCase().trim();

      // log cache misses so they can be collected and added to the cache
      sf_logger.debug("cache miss [{}]", cacheKey + "\t" + nucleotide);
      return nucleotide;
    } catch (Exception e) {
      throw new RuntimeException("Error when requesting DAS data", e);
    }
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
  private String[] getVcfLineFromDefinition(DefinitionFile definitionFile, VariantLocus variantLocus,
      String genomeBuild) {

    //Get star one for Ref column of vcf. Use DAS in future
    NamedAllele referenceAllele = definitionFile.getNamedAlleles().get(0);
    String allele = "";
    // for the exceptions/exemptions
    if (referenceAllele.getAllele(variantLocus) != null) {
      allele = referenceAllele.getAllele(variantLocus);
    }
    int position = variantLocus.getVcfPosition();
    String idField = variantLocus.getRsid();
    if (idField == null) {
      idField = ".";
    }
    String chr = definitionFile.getChromosome();
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
      starAlleles.add(definitionFile.getGeneSymbol() + ":" + namedAllele.getName().replace(" ","") +
          "[" + l.size() + "]" + "is" + namedAllele.getAllele(variantLocus));
    });
    if (alts.size() == 0 ) {
      alts.add(".");
    }
    String refField = expandRepeats(allele);
    // Simplest possible del/ins parsing, presuming the only a single ins or del string, or the correct one being first
    if (alts.size()>0) {
      if (alts.stream().anyMatch(a -> a.contains("ins"))) {
        String nucleotide = getDAS(chr, Integer.toString(position), genomeBuild);
        for (int i = 0; i < alts.size(); i++) {
          alts.set(i, alts.get(i).replace("ins", nucleotide));
        }
        refField = nucleotide;
      }
      if (alts.get(0).contains("del") && !alts.get(0).contains("delGene")) {
        String nucleotide = getDAS(chr, Integer.toString(position), genomeBuild);
        refField = nucleotide + alts.get(0).replace("del", "");
        alts.set(0, nucleotide);
      }
    }
    SortedSet<String> expandedAlts = expandIupacAlts(alts);
    expandedAlts.remove(refField);

    String altField = String.join(",", expandedAlts);
    String infoField;
    if (starAlleles.size() > 0) {
      infoField = "PX=" + String.join(",", starAlleles);
    } else {
      infoField = "POI";
    }

    return new String[]{
        chr,
        Integer.toString(position),
        idField,
        refField,
        altField,
        ".",
        "PASS",
        infoField,
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
        if (!alt.equals("delGene")) {
          bases.add(alt);
        }
      }
    }

    return bases;
  }
}
