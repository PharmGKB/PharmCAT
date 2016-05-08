package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
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
import org.pharmgkb.pharmcat.definition.model.VariantType;


/**
 * Parses {@link DefinitionFile}.
 *
 * @author Ryan Whaley
 */
public class CuratedDefinitionParser {
  public static final String ASSEMBLY = AssemblyMap.GRCH38;
  private static final int sf_minLineCount = 7;
  private static final Splitter sf_tsvSplitter = Splitter.on("\t").trimResults();
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
  private static final Pattern sf_basePattern = Pattern.compile("^(del[ATCG]*)|(ins[ATCG]*)|([ATCGMRWSYKVHDBN]+)|" +
      "([ACTG]+\\([ACTG]+\\)\\d+[ACTG]+)|" +
      "$");
  private static final Pattern sf_hgvsPosition = Pattern.compile("^[cgp]\\.(\\d+).*$");
  private SimpleDateFormat m_dateFormat = new SimpleDateFormat("MM/dd/yy");

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

    parseGeneLine(sf_tsvSplitter.splitToList(lines.get(LINE_GENE)));
    // must parse chromosome first to initialize variants
    parseChromosomeLine(sf_tsvSplitter.splitToList(lines.get(LINE_CHROMOSOME)));
    parseNamingLine(sf_tsvSplitter.splitToList(lines.get(LINE_NAMING)));
    parseProteinLine(sf_tsvSplitter.splitToList(lines.get(LINE_PROTEIN)));
    parseGeneSeqLine(sf_tsvSplitter.splitToList(lines.get(LINE_GENESEQ)));
    parseRSIDLine(sf_tsvSplitter.splitToList(lines.get(LINE_RSID)));
    parsePopLine(sf_tsvSplitter.splitToList(lines.get(LINE_POPS)));
    parseVariantLines(lines);

    // TODO(markwoon): this is temporary, to maintain compatibility with current behavior
    if (m_errors.size() > 0) {
      throw new ParseException(Joiner.on("\n").join(m_errors));
    }

    // find indels and repeats
    for (int x = 0; x < m_definitionFile.getVariants().length; x += 1) {
      VariantLocus locus = m_definitionFile.getVariants()[x];
      NamedAllele refHap = null;
      for (NamedAllele hap : m_definitionFile.getNamedAlleles()) {
        if (refHap == null) {
          refHap = hap;
        }
        String allele = hap.getAlleles()[x];
        if (allele != null) {
          if (allele.contains("del")) {
            if (hap == refHap) {
              locus.setType(VariantType.INS);
              break;
            } else {
              locus.setType(VariantType.DEL);
            }
          }
          if (allele.contains("\\(")) {
            locus.setType(VariantType.REPEAT);
          }
        }
      }
    }

    // finalize NamedAlleles
    for (NamedAllele namedAllele : m_definitionFile.getNamedAlleles()) {
      namedAllele.finalize(m_definitionFile.getVariants());
    }

    return m_definitionFile;
  }

  private void parseGeneLine(List<String> fields) {
    Preconditions.checkNotNull(fields);
    Preconditions.checkArgument(fields.size() >= 2);

    Matcher m = sf_geneFieldPattern.matcher(fields.get(0));
    if (!m.matches()) {
      m_errors.add("Gene field not in expected format.  Expecting 'GENE:geneSymbol', got '" + fields.get(0) + "'");
    } else {
      m_definitionFile.setGeneSymbol(m.group(1));
    }

    try {
      Date date = m_dateFormat.parse(fields.get(1));
      m_definitionFile.setModificationDate(date);
    } catch (java.text.ParseException e) {
      m_errors.add("Couldn't parse date. Expecting " + m_dateFormat.toPattern() + ", got '" + fields.get(1) + "'");
    }
  }

  private void parseChromosomeLine(List<String> fields) {
    String title = fields.get(1);
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

    VariantLocus[] variants = new VariantLocus[fields.size() - 2];
    for (int i = 2; i < fields.size(); i++) {
      String val = StringUtils.stripToNull(fields.get(i));
      if (val != null) {
        variants[i - 2] = parseVariantLocus(i, val);
      }
    }
    m_definitionFile.setVariants(variants);
  }

  private VariantLocus parseVariantLocus(int col, String chrHgvsName) {

    int position = -1;
    for (String hgvs : sf_hgvsSplitter.splitToList(chrHgvsName)) {
      Matcher m = sf_hgvsPosition.matcher(hgvs);
      if (m.matches()) {
        int p =  Integer.parseInt(m.group(1), 10);
        if (position == -1) {
          position = p;
        } else if (position != p) {
          m_errors.add("Cannot handle HGVS name with multiple positions in " + prettyLineCol(LINE_CHROMOSOME, col) +
              ": '" + hgvs + "'");
        }
      } else {
        m_errors.add("Invalid HGVS format in " + prettyLineCol(LINE_CHROMOSOME, col) + ": '" + hgvs + "'");
      }

    }

    return new VariantLocus(position, chrHgvsName);
  }

  private void parseNamingLine(List<String> fields) {
    if (fields.size() <= 2) {
      return;
    }

    for (int i = 2; i < fields.size(); i++) {
      String val = StringUtils.stripToNull(fields.get(i));
      if (val != null) {
        m_definitionFile.getVariants()[i - 2].setResourceNote(val);
      }
    }
  }

  private void parseProteinLine(List<String> fields) {
    String title = fields.get(1);
    if (StringUtils.isBlank(title)) {
      m_errors.add("No protein position description specified in " + prettyLineCol(LINE_PROTEIN, 1));
      return;
    }

    Matcher m = sf_refSeqPattern.matcher(title);
    if (!m.matches()) {
      m_errors.add("Expected protein RefSeq ID in " + prettyLineCol(LINE_PROTEIN, 1));
      return;
    }
    m_definitionFile.setRefSeqProtein(m.group(1));

    for (int i = 2; i < fields.size(); i++) {
      String val = StringUtils.stripToNull(fields.get(i));
      if (val != null) {
        m_definitionFile.getVariants()[i - 2].setProteinNote(fields.get(i));
      }
    }
  }

  private void parseGeneSeqLine(List<String> fields) {

    String title = fields.get(1);
    if (StringUtils.isBlank(title)) {
      m_errors.add("No gene position description specified in " + prettyLineCol(LINE_GENESEQ, 1));
      return;
    }

    Matcher m = sf_refSeqPattern.matcher(title);
    if (!m.matches()) {
      m_errors.add("Expected gene RefSeq ID in " + prettyLineCol(LINE_GENESEQ, 1));
      return;
    }
    m_definitionFile.setRefSeqGene(m.group(1));
    if (title.contains(" forward ")) {
      m_definitionFile.setOrientation("forward");
    } else if (title.contains(" reverse ")) {
      m_definitionFile.setOrientation("reverse");
    } else {
      m_warnings.add("Missing orientation in " + prettyLineCol(LINE_GENESEQ, 1));
    }

    for (int i = 2; i < fields.size(); i++) {
      String val = StringUtils.stripToNull(fields.get(i));
      if (val != null) {
        m_definitionFile.getVariants()[i - 2].setGeneHgvsName(fields.get(i));
      }
    }
  }

  private void parseRSIDLine(List<String> fields) {
    for (int i = 2; i < fields.size(); i++) {
      String val = StringUtils.stripToNull(fields.get(i));
      if (val != null) {
        m_definitionFile.getVariants()[i - 2].setRsid(val);
      }
    }
  }

  private void parsePopLine(List<String> fields) {

    boolean hasError = false;
    if (!fields.get(0).equals("Allele")) {
      m_errors.add("Expected the title 'Allele' in " + prettyLineCol(LINE_POPS, 0));
      hasError = true;
    }
    if (!fields.get(1).equals("Allele Functional Status")) {
      m_errors.add("Expected the title 'Allele Functional Status' in " + prettyLineCol(LINE_POPS, 1));
      hasError = true;
    }
    if (hasError) {
      return;
    }
    if (fields.size() <= 2) {
      m_warnings.add("No populations specified");
      return;
    }

    for (int i = 2; i < fields.size(); i++) {
      String val = StringUtils.stripToNull(fields.get(i));
      if (val != null) {
        Matcher m = sf_populationTitle.matcher(val);
        if (m.matches()) {
          m_definitionFile.addPopulation(val);
          m_populationPositionMap.put(i, val);
        }
      }
    }
  }


  private void parseVariantLines(List<String> lines) {
    boolean isVariantLine = false;
    boolean isNoteLine = false;

    for (int lineNum = 0; lineNum < lines.size(); lineNum += 1) {
      String line = lines.get(lineNum);

      if (isNoteLine) {
        m_definitionFile.addNote(StringUtils.strip(line));

      } else if (line.toLowerCase().startsWith("notes:")) {
        isVariantLine = false;
        isNoteLine = true;

      } else if (isVariantLine) {
        List<String> fields = sf_tsvSplitter.splitToList(line);
        if (fields.size() > 2) {

          String[] variantAlleles = new String[m_definitionFile.getVariants().length];
          for (int i = 2; i < fields.size() && i < variantAlleles.length + 2; i++) {
            String val = StringUtils.stripToNull(fields.get(i));
            if (val != null) {
              if (!sf_basePattern.matcher(val).matches()) {
                m_errors.add("Invalid bases (" + val + ") in " + prettyLineCol(lineNum, i));
              }
              variantAlleles[i - 2] = val;
            }
          }

          String name = StringUtils.stripToNull(fields.get(0));
          String id = m_haplotypeIdMap.get(m_definitionFile.getGeneSymbol(), name);
          if (id == null) {
            m_errors.add("Allele has no ID: " + m_definitionFile.getGeneSymbol() + " " + name);
            return;
          }

          NamedAllele namedAllele = new NamedAllele(id, name, variantAlleles);
          namedAllele.setFunction(fields.get(1));

          Map<String,String> popFreqMap = new HashMap<>();
          m_populationPositionMap.keySet().stream()
              .filter(i -> i < fields.size())
              .forEach(i -> popFreqMap.put(m_populationPositionMap.get(i), fields.get(i)));
          namedAllele.setPopFreqMap(popFreqMap);

          m_definitionFile.addNamedAllele(namedAllele);
        }

      } else if (line.startsWith("Allele")) {
        isVariantLine = true;
      }
    }
  }


  /**
   * Converts a column number to an Excel column name.
   *
   * @param number column number (where 0 = A)
   */
  private static String columnNumberToName(int number) {
    if (number == 0) {
      return "A";
    }
    StringBuilder sb = new StringBuilder();
    while (number > 0) {
      sb.append((char)('A' + (number % 26)));
      number = (number / 26) - 1;
    }
    return sb.reverse().toString();
  }

  /**
   * Formats (0-based) line/column numbers for display.
   *
   * @param line line number
   * @param col column number
   */
  private static String prettyLineCol(int line, int col) {
    return "line " + (line + 1) + ", column " + columnNumberToName(col);
  }
}
