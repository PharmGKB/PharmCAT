package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.apache.commons.lang3.StringUtils;
import org.apache.http.HttpEntity;
import org.apache.http.HttpHeaders;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
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
public class ExtractPositions implements AutoCloseable {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final String sf_refCacheFilename = "positions_reference.tsv";
  private static final Set<String> sf_excludedGenes = ImmutableSet.of("CYP2D6");
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
      "##INFO=<ID=PX,Number=.,Type=String,Description=\"Gene:Allele_Name[n of defining positions]allele\">\n" +
      "##INFO=<ID=POI,Number=0,Type=Flag,Description=\"Position of Interest but not part of an allele definition\">\n" +
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
      "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPharmCAT\n";
  private static final Gson sf_gson = new GsonBuilder()
      .serializeNulls()
      .disableHtmlEscaping()
      .excludeFieldsWithoutExposeAnnotation()
      .setPrettyPrinting().create();

  private final Path m_outputVcf;
  private final DefinitionReader m_definitionReader;
  private final SortedSet<String> m_genes = new TreeSet<>();
  private final CloseableHttpClient m_httpClient = HttpClients.createDefault();
  private final SortedMap<String,String> m_refCache = new TreeMap<>();
  private boolean m_foundNewPosition;


  /**
   * Default constructor.
   *
   * @param outputVcf file path to write output VCF to
   * @param allGenes true to all genes; otherwise will exclude genes for which guidelines are not supported
   */
  public ExtractPositions(Path outputVcf, boolean allGenes) {
    m_outputVcf = outputVcf;
    try {
      // load the allele definition files
      m_definitionReader = new DefinitionReader();
      m_definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
      m_genes.addAll(m_definitionReader.getGenes());
      if (!allGenes) {
        m_genes.removeAll(sf_excludedGenes);
      }
      if (m_genes.size() == 0) {
        throw new RuntimeException("Did not find any allele definitions at " + DataManager.DEFAULT_DEFINITION_DIR);
      }

      // load positions data
      try (BufferedReader reader = Files.newBufferedReader(PathUtils.getPathToResource(getClass(), sf_refCacheFilename))) {
        String line;
        List<String> errors = new ArrayList<>();
        while ((line = StringUtils.stripToNull(reader.readLine())) != null) {
          String[] fields = line.split("\t");
          String previousValue = m_refCache.put(fields[0], fields[1]);
          if (previousValue != null) {
            errors.add("Duplicate key found: " + fields[0]);
          }
        }
        if (errors.size() > 0) {
          throw new RuntimeException("Found problems in the reference cache\n" + String.join("\n", errors));
        }
      }
    } catch (IOException ex) {
      throw new RuntimeException("Error parsing data", ex);
    }
  }

  @Override
  public void close() throws Exception {
    try {
      m_httpClient.close();
    } catch (Exception ex) {
      sf_logger.warn("Error closing http client", ex);
    }

    if (m_foundNewPosition) {
      Path cacheFile = m_outputVcf.getParent().resolve(sf_refCacheFilename);
      System.out.println("\nNew positions found, saving cache file to " + cacheFile);
      System.out.println("Don't forget to save this file!");
      try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(cacheFile))) {
        for (String key : m_refCache.keySet()) {
          writer.println(key + "\t" + m_refCache.get(key));
        }
      } catch (Exception ex) {
        throw new IOException("Error writing new " + sf_refCacheFilename, ex);
      }
    }
  }


  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("o", "output-file", "output vcf file", true, "o")
        .addOption("a", "all-genes", "include all genes");

    cliHelper.execute(args, cli -> {
      try {
        try (ExtractPositions extractPositions =
                 new ExtractPositions(cli.getValidFile("o", false), cli.hasOption("a"))) {
          extractPositions.run();
        }
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
      try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(m_outputVcf))) {
        System.out.println("Writing to " + m_outputVcf);
        writer.print(String.format(sf_fileHeader, DateTimeFormatter.ISO_OFFSET_DATE_TIME.format(ZonedDateTime.now())));

        Map<Integer, Map<Integer, String[]>> chrMap = getPositions();
        for (Integer chr : chrMap.keySet()) {
          Map<Integer, String[]> posMap = chrMap.get(chr);
          for (Integer pos : posMap.keySet()) {
            writer.println(String.join("\t", posMap.get(pos)));
          }
        }
      }
      System.out.println("Done.");
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  /**
   * Build up VCF positions.
   */
  Map<Integer, Map<Integer, String[]>> getPositions() {
    Map<Integer, Map<Integer, String[]>> chrMap = new TreeMap<>();
    for (String gene : m_genes) {
      DefinitionFile definitionFile = m_definitionReader.getDefinitionFile(gene);
      int positionCount = 0;
      // convert bXX format to hgXX
      String build = "hg" + definitionFile.getGenomeBuild().substring(1);
      for (VariantLocus variantLocus : definitionFile.getVariants()) {
        positionCount++;
        addVcfLine(chrMap, getVcfLineFromDefinition(definitionFile, variantLocus, build));
      }
      DefinitionExemption exemption = m_definitionReader.getExemption(gene);
      if (exemption != null) {
        for (VariantLocus variant : exemption.getExtraPositions()) {
          positionCount++;
          addVcfLine(chrMap, getVcfLineFromDefinition(definitionFile, variant, build));
        }
      }
      sf_logger.info("queried {} positions for {}", positionCount, gene);
    }
    return chrMap;
  }

  private void addVcfLine(Map<Integer, Map<Integer, String[]>> chrMap, String[] vcfFields) {
    int chr = Integer.parseInt(vcfFields[0].replace("chr", ""));
    int position = Integer.parseInt(vcfFields[1]);
    Map<Integer, String[]> posMap = chrMap.computeIfAbsent(chr, k -> new TreeMap<>());
    if (posMap.put(position, vcfFields) != null) {
      throw new IllegalStateException("Multiple entries for " + chr + ":" + position);
    }
  }


  /**
   * Helper method to convert repeats into standard format.
   */
  private static String expandRepeats(String allele) {
    String finalAllele = allele;
    if (allele.contains("(")) {
      int bracketStart = 0;
      int bracketEnd = 0;
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
   * Gets the reference allele for the specified position.
   */
  String getReferenceAllele(String assembly, String chr, long position) {
    String cacheKey = assembly + ":" + chr + ":" + position;
    if (m_refCache.containsKey(cacheKey)) {
      return m_refCache.get(cacheKey);
    }

    try {
      String url = "https://api.genome.ucsc.edu/getData/sequence?genome=" + assembly + ";chrom=" + chr + ";start=" +
          (position - 1) + ";end=" + position;
      sf_logger.debug("Getting position: {}", url);
      HttpGet httpGet = new HttpGet(url);
      httpGet.setHeader(HttpHeaders.ACCEPT, "application/xml");
      try (CloseableHttpResponse response = m_httpClient.execute(httpGet)) {
        HttpEntity entity = response.getEntity();
        //noinspection unchecked
        Map<String, String> rez = sf_gson.fromJson(new InputStreamReader(entity.getContent()), Map.class);
        String nucleotide = StringUtils.stripToNull(rez.get("dna"));
        if (nucleotide == null) {
          throw new IOException("UCSC responded with no DNA for " + cacheKey);
        }
        nucleotide = nucleotide.toUpperCase();
        m_refCache.put(cacheKey, nucleotide);
        m_foundNewPosition = true;
        return nucleotide;
      }
    } catch (Exception ex) {
      throw new RuntimeException("Error when requesting DAS data", ex);
    }
  }

  /**
   * Helper method to convert repeats into standard format.
   */
  private String[] getVcfLineFromDefinition(DefinitionFile definitionFile, VariantLocus variantLocus,
      String genomeBuild) {

    long position = variantLocus.getVcfPosition();
    String chr = definitionFile.getChromosome();
    // get alts from the namedAlleles
    List<String> alts = new ArrayList<>();

    String dasRef = getReferenceAllele(genomeBuild, chr, position);
    NamedAllele refNamedAllele = definitionFile.getNamedAlleles().first();
    String referenceAllele;
    // reference named allele may not have allele if the position is an extra or exemption position
    // not used in the allele definition
    if (refNamedAllele.getAllele(variantLocus) != null) {
      referenceAllele = refNamedAllele.getAllele(variantLocus);
      if (referenceAllele.length() == 1 && Iupac.lookup(referenceAllele).isAmbiguity()) {
        alts.addAll(Iupac.lookup(referenceAllele).getBases());
        if (alts.remove(dasRef)) {
          referenceAllele = dasRef;
        } else {
          // this will be undone in dasRef check below, but do it anyway so we don't have to duplicate warning message
          referenceAllele = alts.get(0);
          alts.remove(referenceAllele);
        }
      }
      if (!dasRef.equals(referenceAllele) && variantLocus.getType() == VariantType.SNP && referenceAllele.length() == 1) {
        System.err.println(definitionFile.getGeneSymbol() + " " + refNamedAllele.getName() +
            " has ref mismatch at position " + position +
            (variantLocus.getRsid() != null ? " (" + variantLocus.getRsid() + ")" : "") + " - using " + dasRef +
            " as ref as opposed to " + referenceAllele);
        alts.add(referenceAllele);
        referenceAllele = dasRef;
      }
    } else {
      referenceAllele = dasRef;
    }

    List<String> pxInfo = new ArrayList<>();
    for (NamedAllele namedAllele : definitionFile.getNamedAlleles()) {
      String allele = namedAllele.getAllele(variantLocus);
      if (allele == null) {
        continue;
      }
      List<String> alleles;;
      if (allele.length() == 1) {
        alleles = Iupac.lookup(allele).getBases();
      } else {
        alleles = Lists.newArrayList(expandRepeats(allele));
      }
      for (String a : alleles) {
        if (!alts.contains(a) && !referenceAllele.equals(a)) {
          alts.add(a);
        }
      }
      long numAlleles = Arrays.stream(namedAllele.getAlleles())
          .filter(Objects::nonNull)
          .count();
      pxInfo.add(definitionFile.getGeneSymbol() + ":" + namedAllele.getName().replace(" ","") +
          "[" + numAlleles + "]" + "is" + namedAllele.getAllele(variantLocus));
    }
    if (alts.size() == 0 ) {
      alts.add(".");
    }

    String refField = expandRepeats(referenceAllele);
    // simplest possible del/ins parsing, presuming the only a single ins or del string, or the correct one being first
    if (alts.size() > 0) {
      if (alts.stream().anyMatch(a -> a.contains("ins"))) {
        String nucleotide = getReferenceAllele(genomeBuild, chr, position);
        for (int i = 0; i < alts.size(); i++) {
          alts.set(i, alts.get(i).replace("ins", nucleotide));
        }
        refField = nucleotide;
      }
      if (alts.get(0).contains("del") && !alts.get(0).contains("delGene")) {
        String nucleotide = getReferenceAllele(genomeBuild, chr, position);
        refField = nucleotide + alts.get(0).replace("del", "");
        alts.set(0, nucleotide);
      }
      // sanity check
      if (alts.stream().anyMatch(a -> a.contains("del"))) {
        throw new IllegalStateException("An ALT allele contains 'del': " + String.join("/", alts));
      }
    }

    String idField = variantLocus.getRsid();
    if (idField == null) {
      idField = ".";
    }

    String infoField;
    if (pxInfo.size() > 0) {
      infoField = "PX=" + String.join(",", pxInfo);
    } else {
      infoField = "POI";
    }

    return new String[]{
        chr,
        Long.toString(position),
        idField,
        refField,
        String.join(",", new TreeSet<>(alts)),
        ".",
        "PASS",
        infoField,
        "GT",
        "0/0"
    };
  }

  /**
   * Gets the genes applicable to this process which excludes genes we are not currently supporting for allele matching.
   *
   * @return a Set of gene symbols
   */
  Set<String> getGenes() {
    return m_genes;
  }
}
