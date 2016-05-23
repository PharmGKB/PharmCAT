package org.pharmgkb.pharmcat.reporter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Set;
import javax.annotation.Nonnull;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.gson.Gson;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.reporter.model.CPICException;
import org.pharmgkb.pharmcat.reporter.model.CPICExceptionList;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;


/**
 *
 * Core class for matching called genes and haplotypes to the exceptions list
 *
 * @author greytwist
 */

public class ExceptionMatcher {
  private Multimap<String, CPICException> m_exceptionMap;

  public ExceptionMatcher() throws Exception {
    try {
      File exceptionFile = new File(getClass().getResource("exceptions.json").toURI());
      m_exceptionMap = loadExceptions(exceptionFile);
    } catch (IOException|URISyntaxException e) {
      throw new Exception("Not able to parse excepiton list", e);
    }
  }

  /**
   * Adds any matching exceptions to the given {@link GeneReport}
   */
  public void addExceptions(@Nonnull GeneReport geneReport) {
    String geneSymbol = geneReport.getGene();

    if (m_exceptionMap.containsKey(geneSymbol) ){
      // add any known gene exceptions
      m_exceptionMap.get(geneSymbol).stream()
          .filter(exception -> test(geneReport, exception.getMatches()))
          .forEach(geneReport::addException);
    }
  }

  /**
   * FIXME: here I make the assumption of three different rule types
   * this should be drivin in config or some other place that makes handleing
   * matching easy without having to edit code
   */
  public boolean test(GeneReport gene, String matchLogic) {
    boolean finalTest = false;
    boolean geneTest = false;
    boolean hapTest = false;
    boolean dipTest = false;

    String[] parts = matchLogic.replaceAll("\\s", "").split(",");
    for (String part : parts) {
      String[] type = part.split("=");
      if (type[0].equals("gene")) {
        geneTest = true; // just assume this is true for now, should double check that gene strings match
      } else if (type[0].equals("hap")) {
        if( type[1].replaceAll("'", "").isEmpty()){
          hapTest = true;
        } else {
          String cleanHap = type[1].replaceAll("'", "");
          hapTest = testForHapCall(cleanHap, gene.getDips() );
        }

      } else if( type[0].equals("dip") ) {
        if( type[1].replaceAll("'", "").isEmpty()){
          dipTest = true;
        } else {
          String cleanDip = type[1].replaceAll("'", "");
          dipTest = testForDipCall( cleanDip, gene.getDips() );
        }
      }
    }
    if( geneTest && hapTest && dipTest ){
      finalTest = true;
    }

    return finalTest;
  }

  private boolean testForDipCall(String cleanDip, Set<String> dips) {
    boolean test = false;
    for( String dip : dips) {
      if( HaplotypeNameComparator.getComparator().compare(cleanDip, dip) == 0 ){
        test = true;
        break;
      }
    }
    return test;
  }

  private boolean testForHapCall(String cleanHap, Set<String> dips) {
    boolean test = false;
    for( String dip : dips) {
      String[] parts = dip.split("/");
      if( parts[0].equals(cleanHap) || parts[1].equals(cleanHap)){
        test = true;
        break;
      }
    }
    return test;
  }

  private Multimap<String, CPICException> loadExceptions(File exceptions)throws IOException {
    Gson gson = new Gson();
    Multimap<String, CPICException> matcher = HashMultimap.create();

    try (BufferedReader br = new BufferedReader(new FileReader(exceptions))) {
      CPICExceptionList exceptionList = gson.fromJson(br, CPICExceptionList.class);

      for (CPICException rule : exceptionList.getRules()) {
        matcher.put(rule.getGene(), rule);
      }
    }
    return matcher;
  }
}
