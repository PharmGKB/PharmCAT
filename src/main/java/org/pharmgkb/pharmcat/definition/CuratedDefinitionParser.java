package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.google.common.base.Joiner;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;


/**
 * Parses {@link DefinitionFile}.
 *
 * @author Ryan Whaley
 */
public class CuratedDefinitionParser {
  public static final String ASSEMBLY = AssemblyMap.GRCH38;
  private static final int sf_minLineCount = 7;
  private static final String sf_separator = "\t";
  private static final Splitter sf_hgvsSplitter = Splitter.on(";").trimResults();

  private static final int LINE_GENE = 0;
  private static final int LINE_NAMING = 1;
  private static final int LINE_PROTEIN = 2;
  private static final int LINE_CHROMOSOME = 3;
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
  private DefinitionFile m_definitionFile;
  private Map<Integer,String> m_populationPositionMap = new HashMap<>();
  private List<String> m_warnings = new ArrayList<>();
  private List<String> m_errors = new ArrayList<>();
  private HaplotypeIdMap m_haplotypeIdMap;


  public CuratedDefinitionParser(Path filePath) throws IOException {
    Preconditions.checkNotNull(filePath);
    Preconditions.checkArgument(filePath.toFile().exists());
    Preconditions.checkArgument(filePath.toFile().isFile());
    m_filePath = filePath;
    m_definitionFile = new DefinitionFile();
    m_haplotypeIdMap = new HaplotypeIdMap();
  }

  public List<String> getWarnings() {
    return m_warnings;
  }

  public DefinitionFile parse() {

    List<String> lines;
    try {
      lines = Files.readAllLines(m_filePath);
    } catch (IOException e) {
      throw new ParseException("Couldn't read lines in "+m_filePath, e);
    }
    Preconditions.checkState(lines.size() >= sf_minLineCount, "Not enough lines in the translation file");

    parseGeneLine(lines.get(LINE_GENE).split(sf_separator));
    parseChromosomeLine(lines.get(LINE_CHROMOSOME).split(sf_separator));
    parseNamingLine(lines.get(LINE_NAMING).split(sf_separator));
    parseProteinLine(lines.get(LINE_PROTEIN).split(sf_separator));
    parseGeneSeqLine(lines.get(LINE_GENESEQ).split(sf_separator));
    parseRSIDLine(lines.get(LINE_RSID).split(sf_separator));
    parsePopLine(lines.get(LINE_POPS).split(sf_separator));
    parseVariantLines(lines);

    // TODO(markwoon): this is temporary, to maintain compatibility with current behavior
    if (m_errors.size() > 0) {
      throw new ParseException(Joiner.on("\n").join(m_errors));
    }

    for (int x = 0; x < m_definitionFile.getVariants().length; x += 1) {
      VariantLocus locus = m_definitionFile.getVariants()[x];
      for (NamedAllele namedAllele : m_definitionFile.getNamedAlleles()) {
        String allele = namedAllele.getAlleles()[x];
        if (allele != null && (allele.contains("ins") || allele.contains("del"))) {
          locus.setInDel(true);
          break;
        }
      }
    }

    return m_definitionFile;
  }

  private void parseGeneLine(String[] fields) {
    Preconditions.checkNotNull(fields);
    Preconditions.checkArgument(fields.length >= 2);

    Matcher m = sf_geneFieldPattern.matcher(fields[0]);
    if (!m.matches()) {
      m_errors.add("Gene field not in expected format.  Expecting 'GENE:geneSymbol', got '" + fields[0] + "'");
    } else {
      m_definitionFile.setGeneSymbol(m.group(1));
    }

    try {
      Date date = sf_dateFormat.parse(fields[1]);
      m_definitionFile.setModificationDate(date);
    } catch (java.text.ParseException e) {
      m_errors.add("Couldn't parse date. Expecting " + sf_dateFormat.toPattern() + ", got '" + fields[1] + "'");
    }
  }

  private void parseChromosomeLine(String[] fields) {
    String title = fields[1];
    Preconditions.checkState(StringUtils.isNotBlank(title), "No chromosomal position description specified");

    Matcher m = sf_refSeqPattern.matcher(title);
    if (!m.matches()) {
      m_errors.add("Missing chromosome RefSeq ID in line " + LINE_CHROMOSOME);

    } else {
      m_definitionFile.setRefSeqChromosome(m.group(1));

      AssemblyMap assemblyMap;
      try {
        assemblyMap = new AssemblyMap();
      } catch (IOException e) {
        throw new RuntimeException("Couldn't make assembly to build map", e);
      }
      String build = assemblyMap.get(m_definitionFile.getRefSeqChromosome());
      if (build == null) {
        m_errors.add(m_definitionFile.getRefSeqChromosome() + " not associated with any build");
      } else if (!build.equals(ASSEMBLY)) {
        m_errors.add("Chromosome identifier not on GRCh " + ASSEMBLY + ": " + m_definitionFile.getRefSeqChromosome());
      } else {
        m_definitionFile.setGenomeBuild(build);
      }

      int chrNum = Integer.parseInt(m.group(2), 10); // a leading 0 sometimes indicates octal, but we know this is always base 10
      if (!(chrNum >= 1 && chrNum <= 24)) {
        throw new RuntimeException("Unknown or unsupported chromosome number " + chrNum + " on chromosomal line " + LINE_CHROMOSOME);
      }
      String chromosomeName;
      if (chrNum == 23) {
        chromosomeName = "chrX";
      } else if (chrNum == 24) {
        chromosomeName = "chrY";
      } else {
        chromosomeName = "chr" + chrNum;
      }
      m_definitionFile.setChromosome(chromosomeName);
    }

    m_definitionFile.setVariants(new VariantLocus[fields.length-2]);
    for (int i=2; i<fields.length; i++) {
      if (StringUtils.isNotBlank(fields[i])) {
        VariantLocus variantLocus = parseVariantLocus(i, fields[i]);
        m_definitionFile.getVariants()[i-2] = variantLocus;
      }
    }
  }

  private void parseNamingLine(String[] fields) {
    if (fields.length <= 2) {
      return;
    }

    for (int i=2; i<fields.length; i++) {
      m_definitionFile.getVariants()[i-2].setResourceNote(fields[i]);
    }
  }

  private void parseProteinLine(String[] fields) {
    String title = fields[1];
    Preconditions.checkState(StringUtils.isNotBlank(title), "No protein position description specified");

    Matcher m = sf_refSeqPattern.matcher(title);
    if (!m.matches()) {
      throw new RuntimeException("Expected protein RefSeq ID in line "+LINE_PROTEIN);
    }
    m_definitionFile.setRefSeqProtein(m.group(1));

    for (int i=2; i<fields.length; i++) {
      m_definitionFile.getVariants()[i-2].setProteinNote(fields[i]);
    }
  }

  private VariantLocus parseVariantLocus(int col, String text) {
    VariantLocus variantLocus = new VariantLocus();

    sf_hgvsSplitter.splitToList(text).stream()
        .forEach(p -> {
          if (!sf_hgvsPosition.matcher(StringUtils.strip(p)).matches()) {
            m_errors.add("Invalid HGVS position format in column " + columnNumberToName(col) + ": '" + p + "'");
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
      throw new RuntimeException("Expected gene RefSeq ID in line " + LINE_GENESEQ);
    }
    m_definitionFile.setRefSeqGene(m.group(1));
    if (title.contains(" forward ")) {
      m_definitionFile.setOrientation("forward");
    } else if (title.contains(" reverse ")) {
      m_definitionFile.setOrientation("reverse");
    } else {
      m_warnings.add("Missing orientation");
    }

    for (int i=2; i<fields.length; i++) {
      m_definitionFile.getVariants()[i-2].setGenePosition(fields[i]);
    }
  }

  private void parseRSIDLine(String[] fields) {
    for (int i=2; i < fields.length; i++) {
      m_definitionFile.getVariants()[i-2].setRsid(fields[i]);
    }
  }

  private void parsePopLine(String[] fields) {
    Preconditions.checkArgument(fields[0].equals("Allele"),
        "Expected the title 'Allele' in first column of row " + (LINE_POPS+1));
    Preconditions.checkArgument(fields[1].equals("Allele Functional Status"),
        "Expected the title 'Allele Functional Status' in second column of row " + (LINE_POPS+1));
    Preconditions.checkArgument(fields.length>2, "No populations specified");

    for (int i=2; i<fields.length; i++) {
      String pop = fields[i];
      Matcher m = sf_populationTitle.matcher(pop);
      if (m.matches()) {
        m_definitionFile.addPopulation(pop);
        m_populationPositionMap.put(i, pop);
      }
    }
  }


  private void parseVariantLines(List<String> lines) {
    boolean isVariantLine = false;
    boolean isNoteLine = false;
    int lastVariantIndex = 2 + m_definitionFile.getVariants().length;

    for (String line : lines) {
      if (isNoteLine) {
        m_definitionFile.addNote(StringUtils.strip(line));

      } else if (line.toLowerCase().startsWith("notes:")) {
        isVariantLine = false;
        isNoteLine = true;

      } else if (isVariantLine) {
        String[] fields = line.split(sf_separator);
        if (fields.length > 2) {
          String[] variantAlleles = Arrays.copyOfRange(fields, 2, lastVariantIndex);
          for (int i = 0; i < variantAlleles.length; i++) {
            variantAlleles[i] = StringUtils.stripToNull(variantAlleles[i]);
          }
          validateAlleles(fields[0], variantAlleles);

          NamedAllele namedAllele = new NamedAllele();
          namedAllele.setName(fields[0]);
          namedAllele.setFunction(fields[1]);
          namedAllele.setAlleles(variantAlleles);

          Map<String,String> popFreqMap = new HashMap<>();
          m_populationPositionMap.keySet().stream()
              .filter(i -> i < fields.length)
              .forEach(i -> popFreqMap.put(m_populationPositionMap.get(i), fields[i]));
          namedAllele.setPopFreqMap(popFreqMap);

          String id = m_haplotypeIdMap.get(m_definitionFile.getGeneSymbol(), namedAllele.getName());
          if (id == null) {
            throw new ParseException("Allele has no ID: " + m_definitionFile.getGeneSymbol() + " " + namedAllele.getName());
          }
          namedAllele.setId(id);

          m_definitionFile.addNamedAllele(namedAllele);
        }

      } else if (line.startsWith("Allele")) {
        isVariantLine = true;
      }
    }
  }

  private void validateAlleles(String alleleName, String[] alleles) {

    StringBuilder errBuilder = new StringBuilder();
    for (String allele : alleles) {
      if (allele != null && !sf_basePattern.matcher(allele).matches()) {
        if (errBuilder.length() > 0){
          errBuilder.append("; ");
        }
        errBuilder.append(allele);
      }
    }
    if (errBuilder.length() > 0) {
      throw new ParseException("Allele " + alleleName + " has bad base pair values: " + errBuilder.toString());
    }
  }


  private static String columnNumberToName(int number) {
    StringBuilder sb = new StringBuilder();
    while (number > 0) {
      sb.append((char)('A' + (number % 26)));
      number = (number / 26) - 1;
    }
    return sb.reverse().toString();
  }
}
