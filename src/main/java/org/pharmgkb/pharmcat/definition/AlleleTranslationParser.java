package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.definition.model.AlleleTranslation;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;

/**
 * @author Ryan Whaley
 */
public class AlleleTranslationParser {
  private static final int sf_minLineCount = 7;
  private static final String sf_separator = "\t";
  private static final String sf_hgvsSeparator = ";";

  private static final int LINE_GENE = 0;
  private static final int LINE_NAMING = 1;
  private static final int LINE_PROTEIN = 2;
  private static final int LINE_CHROMO = 3;
  private static final int LINE_GENESEQ = 4;
  private static final int LINE_RSID = 5;
  private static final int LINE_POPS = 6;

  private static final Pattern sf_geneFieldPattern = Pattern.compile("^GENE:\\s*(\\w+)$");
  private static final Pattern sf_refSeqPattern = Pattern.compile("^.*(N\\w_(\\d+)\\.\\d+).*$");
  private static final Pattern sf_populationTitle = Pattern.compile("^(.*) Allele Frequency$");
  private static final Pattern sf_basePattern = Pattern.compile("^(del[ATCG]*)|(ins[ATCG]*)|([ATCGMRWSYKVHDBN]*)$");
  private static final Pattern sf_hgvsPosition = Pattern.compile("^[cgp]\\.(\\d+).*$");
  private static final SimpleDateFormat sf_dateFormat = new SimpleDateFormat("MM/dd/yy");

  private Path m_filePath;
  private AlleleTranslation m_alleleTranslation;
  private Map<Integer,String> m_populationPositionMap = new HashMap<>();

  public AlleleTranslationParser(Path filePath) {
    Preconditions.checkNotNull(filePath);
    Preconditions.checkArgument(filePath.toFile().exists());
    Preconditions.checkArgument(filePath.toFile().isFile());
    m_filePath = filePath;
    m_alleleTranslation = new AlleleTranslation();
  }

  public AlleleTranslation parse() {
    List<String> lines;
    try {
      lines = Files.readAllLines(m_filePath);
    } catch (IOException e) {
      throw new RuntimeException("Couldn't read lines in "+m_filePath, e);
    }
    Preconditions.checkState(lines.size() >= sf_minLineCount, "Not enough lines in the translation file");

    parseGeneLine(lines.get(LINE_GENE).split(sf_separator));
    parseChromoLine(lines.get(LINE_CHROMO).split(sf_separator));
    parseNamingLine(lines.get(LINE_NAMING).split(sf_separator));
    parseProteinLine(lines.get(LINE_PROTEIN).split(sf_separator));
    parseGeneSeqLine(lines.get(LINE_GENESEQ).split(sf_separator));
    parseRSIDLine(lines.get(LINE_RSID).split(sf_separator));
    parsePopLine(lines.get(LINE_POPS).split(sf_separator));
    parseVariantLines(lines);

    return m_alleleTranslation;
  }

  private void parseGeneLine(String[] fields) {
    Preconditions.checkNotNull(fields);
    Preconditions.checkArgument(fields.length >= 2);

    Matcher m = sf_geneFieldPattern.matcher(fields[0]);
    if (!m.matches()) {
      throw new RuntimeException("Gene field not in expected format: "+fields[0]);
    }

    m_alleleTranslation.setGeneSymbol(m.group(1));

    try {
      Date date = sf_dateFormat.parse(fields[1]);
      m_alleleTranslation.setModificationDate(date);
    } catch (ParseException e) {
      throw new RuntimeException("Couldn't parse date "+fields[1], e);
    }
  }

  private void parseChromoLine(String[] fields) {
    String title = fields[1];
    Preconditions.checkState(StringUtils.isNotBlank(title), "No chromosomal position description specified");

    Matcher m = sf_refSeqPattern.matcher(title);
    if (!m.matches()) {
      throw new RuntimeException("Expected chromsome RefSeq ID in line "+LINE_CHROMO);
    }
    m_alleleTranslation.setRefSeqChromo(m.group(1));

    AssemblyMap assemblyMap;
    try {
      assemblyMap = new AssemblyMap();
    } catch (IOException e) {
      throw new RuntimeException("Couldn't make assembly to build map", e);
    }
    String build = assemblyMap.get(m_alleleTranslation.getRefSeqChromo());
    if (build == null) {
      throw new RuntimeException(m_alleleTranslation.getRefSeqChromo() + " not associated with any build");
    }
    if (!build.equals(AssemblyMap.GRCH38)) {
      throw new RuntimeException("Chromosome identifier not on GRCh38: " + m_alleleTranslation.getRefSeqChromo());
    }
    m_alleleTranslation.setGenomeBuild(build);

    int chrNum = Integer.parseInt(m.group(2), 10); // a leading 0 sometimes indicates octal, but we know this is always base 10
    if (!(chrNum >= 1 && chrNum <= 24)) {
      throw new RuntimeException("Unknown or unsupported chromosome number "+chrNum+" on chromosomal line "+LINE_CHROMO);
    }
    String chromosomeName;
    if (chrNum == 23) {
      chromosomeName = "chrX";
    } else if (chrNum == 24) {
      chromosomeName = "chrY";
    } else {
      chromosomeName = "chr" + chrNum;
    }
    m_alleleTranslation.setChromoName(chromosomeName);

    m_alleleTranslation.setVariants(new VariantLocus[fields.length-2]);
    for (int i=2; i<fields.length; i++) {
      if (StringUtils.isNotBlank(fields[i])) {
        VariantLocus variantLocus = parseVariantLocus(fields[i]);
        m_alleleTranslation.getVariants()[i-2] = variantLocus;
      }
    }
  }

  private void parseNamingLine(String[] fields) {
    if (fields.length <= 2) {
      return;
    }

    for (int i=2; i<fields.length; i++) {
      m_alleleTranslation.getVariants()[i-2].setResourceNote(fields[i]);
    }
  }

  private void parseProteinLine(String[] fields) {
    String title = fields[1];
    Preconditions.checkState(StringUtils.isNotBlank(title), "No protein position description specified");

    Matcher m = sf_refSeqPattern.matcher(title);
    if (!m.matches()) {
      throw new RuntimeException("Expected protein RefSeq ID in line "+LINE_PROTEIN);
    }
    m_alleleTranslation.setRefSeqProtein(m.group(1));

    for (int i=2; i<fields.length; i++) {
      m_alleleTranslation.getVariants()[i-2].setProteinNote(fields[i]);
    }
  }

  private VariantLocus parseVariantLocus(String text) {
    VariantLocus variantLocus = new VariantLocus();

    String[] positions = text.split(sf_hgvsSeparator);
    Arrays.stream(positions).forEach(p -> {
      if (!sf_hgvsPosition.matcher(StringUtils.strip(p)).matches()) {
        throw new RuntimeException("Position format doesn't match: "+p);
      }
    });

    variantLocus.setChrPosition(text);
    return variantLocus;
  }

  private void parseGeneSeqLine(String[] fields) {
    String title = fields[1];
    Preconditions.checkState(StringUtils.isNotBlank(title), "No gene position description specified");

    Matcher m = sf_refSeqPattern.matcher(title);
    if (!m.matches()) {
      throw new RuntimeException("Expected gene RefSeq ID in line "+LINE_GENESEQ);
    }
    m_alleleTranslation.setRefSeqGene(m.group(1));

    for (int i=2; i<fields.length; i++) {
      m_alleleTranslation.getVariants()[i-2].setGenePosition(fields[i]);
    }
  }

  private void parseRSIDLine(String[] fields) {
    for (int i=2; i < fields.length; i++) {
      m_alleleTranslation.getVariants()[i-2].setRsid(fields[i]);
    }
  }

  private void parsePopLine(String[] fields) {
    Preconditions.checkArgument(fields[0].equals("Allele"), "Expected the title 'Allele' in first column of row " + (LINE_POPS+1));
    Preconditions.checkArgument(fields[1].equals("Allele Functional Status"), "Expected the title 'Allele Functional Status' in second column of row " + (LINE_POPS+1));
    Preconditions.checkArgument(fields.length>2, "No populations specified");

    for (int i=2; i<fields.length; i++) {
      String pop = fields[i];
      Matcher m = sf_populationTitle.matcher(pop);
      if (m.matches()) {
        m_alleleTranslation.addPopulation(pop);
        m_populationPositionMap.put(i, pop);
      }
    }
  }

  private void parseVariantLines(List<String> lines) {
    HaplotypeIdMap geneHapMap;
    try {
      geneHapMap = new HaplotypeIdMap();
    } catch (IOException e) {
      throw new RuntimeException("Couldn't load haplotype ID map", e);
    }
    boolean isVariantLine = false;
    boolean isNoteLine = false;
    int lastVariantIndex = 2 + m_alleleTranslation.getVariants().length;

    for (String line : lines) {
      if (isNoteLine) {
        m_alleleTranslation.addNote(StringUtils.strip(line));
      }
      else if (line.toLowerCase().startsWith("notes:")) {
        isVariantLine = false;
        isNoteLine = true;
      }
      else if (isVariantLine) {
        String[] fields = line.split(sf_separator);
        if (fields.length > 2) {
          ArrayList<String> alleles = parseVariantFields(fields);
          Set<String> badAlleles = alleles
              .stream()
              .filter(f -> !sf_basePattern.matcher(f).matches())
              .collect(Collectors.toSet());
          if (badAlleles.size()>0) {
            throw new RuntimeException(fields[0] + " has bad base pair values " + badAlleles.stream().collect(Collectors.joining(";")));
          }

          NamedAllele namedAllele = new NamedAllele();
          namedAllele.setName(fields[0]);
          namedAllele.setFunction(fields[1]);

          String[] variantAlleles = Arrays.copyOfRange(fields, 2, lastVariantIndex);
          for (int i = 0; i < variantAlleles.length; i++) {
            variantAlleles[i] = StringUtils.stripToNull(variantAlleles[i]);
          }
          namedAllele.setAlleles(variantAlleles);

          Map<String,String> popFreqMap = new HashMap<>();
          m_populationPositionMap.keySet().stream()
              .filter(i -> i < fields.length)
              .forEach(i -> popFreqMap.put(m_populationPositionMap.get(i), fields[i]));
          namedAllele.setPopFreqMap(popFreqMap);


          String id = geneHapMap.get(m_alleleTranslation.getGeneSymbol(), namedAllele.getName());
          if (id == null) {
            throw new RuntimeException("Allele has no ID: "+m_alleleTranslation.getGeneSymbol()+" "+namedAllele.getName());
          }
          namedAllele.setId(id);

          m_alleleTranslation.addNamedAllele(namedAllele);
        }
      }
      else if (line.startsWith("Allele")) {
        isVariantLine = true;
      }
    }
  }

  private ArrayList<String> parseVariantFields(String[] fields) {
    ArrayList<String> list = new ArrayList<>();
    for (int i = 2;  i < 2 + m_alleleTranslation.getVariants().length;  i++) {
      if (i < fields.length) {
        list.add(fields[i]);
      } else {
        list.add("");
      }
    }
    return list;
  }
}
