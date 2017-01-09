package org.pharmgkb.pharmcat.reporter;

import java.io.BufferedReader;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.gson.Gson;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.reporter.model.CPICException;
import org.pharmgkb.pharmcat.reporter.model.CPICExceptionList;
import org.pharmgkb.pharmcat.reporter.model.CPICExceptionMatch;
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
      m_exceptionMap = loadExceptions(getClass().getResource("exceptions.json").toURI());
    } catch (IOException|URISyntaxException e) {
      throw new Exception("Not able to parse exception list", e);
    }
  }

  /**
   * Adds any matching exceptions to the given {@link GeneReport}
   * RMW: this needs to be fixed
   */
  public void addExceptions(@Nonnull GeneReport geneReport) {
    String geneSymbol = geneReport.getGene();

    if (m_exceptionMap.containsKey(geneSymbol)) {
      // add any known gene exceptions
      m_exceptionMap.get(geneSymbol).stream()
          .filter(exception -> test(geneReport, exception))
          .forEach(geneReport::addException);
    }
  }

  /**
   * testing the observed allele data against the exception logic
   */
  public boolean test(GeneReport geneReport, CPICException exception) {

    CPICExceptionMatch matchCriteria = exception.getMatches();

    // if gene doesn't match, nothing else matters
    if (!matchCriteria.getGenes().contains(geneReport.getGene())) {
      return false;
    }

    // if we're only matching on Gene then we're done
    if (matchCriteria.getHapsCalled().isEmpty() && matchCriteria.getDips().isEmpty()) {
      return true;
    }

    if (geneReport.getMatchData() != null && !geneReport.getMatchData().getHaplotypes().isEmpty()) {
      List<String> observedHaps = geneReport.getMatchData().getHaplotypes().stream()
          .map(NamedAllele::getName)
          .collect(Collectors.toList());
      return testForHapCall(matchCriteria.getHapsCalled(), observedHaps);
    }

    return geneReport.getMatchData() != null 
        && !geneReport.getDips().isEmpty() 
        && testForDipCall(matchCriteria.getDips(), geneReport.getDips());

  }

  private boolean testForDipCall(Collection<String> cleanDip, Collection<String> dips) {
    return !Collections.disjoint(cleanDip, dips);
  }

  private boolean testForHapCall(Collection<String> testHaplotypes, Collection<String> observedHaplotypes) {
    return !Collections.disjoint(testHaplotypes, observedHaplotypes);
  }

  private Multimap<String, CPICException> loadExceptions(URI exceptionsUri)throws IOException {
    Gson gson = new Gson();
    Multimap<String, CPICException> matcher = HashMultimap.create();

    try (BufferedReader br = Files.newBufferedReader(Paths.get(exceptionsUri))) {
      CPICExceptionList exceptionList = gson.fromJson(br, CPICExceptionList.class);

      for (CPICException rule : exceptionList.getRules()) {
        rule.getMatches().getGenes().forEach(g -> matcher.put(g, rule));
      }
    }
    return matcher;
  }
}
