package org.pharmgkb.pharmcat.reporter;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;


/**
 * Factory class for spitting out {@link Diplotype} objects from known information about the gene
 *
 * This is designed to take information from either the Matcher or from outside calls and produce consistent results
 *
 * @author Ryan Whaley
 */
public class DiplotypeFactory {
  private static final Set<String> PHENOTYPE_ONLY = ImmutableSet.of("HLA-A", "HLA-B");
  private static final String UNASSIGNED_FUNCTION = "Unassigned function";

  private final String f_gene;
  private final String f_referenceAlleleName;
  private final GenePhenotype f_genePhenotype;
  private final Map<String,Haplotype> m_haplotypeCache = new HashMap<>();

  /**
   * public constructor
   *
   * Initialize the factory with all the necessary information to make haplotype and diplotype calls
   *
   * @param gene the gene symbol for the diplotypes this will call
   * @param genePhenotype a {@link GenePhenotype} object that maps haplotypes and functions to phenotypes
   * @param referenceAlleleName the name of the reference allele
   */
  public DiplotypeFactory(String gene, GenePhenotype genePhenotype, @Nullable String referenceAlleleName) {
    f_gene = gene;
    f_genePhenotype = genePhenotype;
    f_referenceAlleleName = referenceAlleleName;
  }

  /**
   * Make diplotype objects based on GeneCall objects that come from the NamedAlleleMatcher
   */
  public List<Diplotype> makeDiplotypes(GeneCall geneCall) {
    Preconditions.checkNotNull(geneCall);
    Preconditions.checkArgument(geneCall.getGene().equals(f_gene));

    // do the regular processing when diplotypes are called
    if (geneCall.getDiplotypes().size() > 0) {
      return geneCall.getDiplotypes().stream().map(this::makeDiplotype).collect(Collectors.toList());
    }

    // if not diplotypes are matched then mark it as unknown
    else {
      return ImmutableList.of(makeUnknownDiplotype());
    }
  }

  /**
   * Make diplotype objects based on {@link OutsideCall} objects
   */
  public List<Diplotype> makeDiplotypes(OutsideCall outsideCall) {
    Preconditions.checkNotNull(outsideCall);
    Preconditions.checkArgument(outsideCall.getGene().equals(f_gene));

    if (outsideCall.getDiplotypes().size() > 0) {
      return makeDiplotypes(outsideCall.getDiplotypes());
    } else {
      return makeDiplotypesFromPhenotype(outsideCall.getPhenotype());
    }
  }

  /**
   * Make Diplotype objects based on string diplotype names like "*1/*80"
   * @param dipNames a collection of diplotype names like "*1/*80"
   * @return a List of full Diplotype objects
   */
  public List<Diplotype> makeDiplotypes(@Nullable Collection<String> dipNames) {
    if (dipNames == null) {
      return new ArrayList<>();
    }
    return dipNames.stream().map(this::makeDiplotype).collect(Collectors.toList());
  }

  public List<Diplotype> makeDiplotypesFromPhenotype(String phenotype) {
    Diplotype diplotype = new Diplotype(f_gene, phenotype);

    return ImmutableList.of(diplotype);
  }

  /**
   * Make a Diplotype object based on a DiplotypeMatch, this will also populate function information from the match
   */
  private Diplotype makeDiplotype(DiplotypeMatch diplotypeMatch) {
    HaplotypeMatch h1 = diplotypeMatch.getHaplotype1();
    HaplotypeMatch h2 = diplotypeMatch.getHaplotype2();

    Diplotype diplotype = new Diplotype(f_gene, makeHaplotype(h1), makeHaplotype(h2));
    fillDiplotype(diplotype);

    return diplotype;
  }

  /**
   * Make a Diplotype based on a string in the form "*1/*20"
   * @param diplotypeText a string in the form "*1/*20"
   * @return a Diplotype containing the specified haplotypes
   */
  private Diplotype makeDiplotype(String diplotypeText) {
    Preconditions.checkArgument(StringUtils.isNotBlank(diplotypeText));

    Diplotype diplotype;
    if (diplotypeText.contains("/")) {
      if (GeneReport.isSinglePloidy(f_gene) && !GeneReport.isXChromo(f_gene)) {
        throw new RuntimeException("Cannot specify two genotypes [" + diplotypeText + "] for single chromosome gene " + f_gene);
      }

      String[] alleles = diplotypeText.split("/");
      diplotype = new Diplotype(f_gene, makeHaplotype(alleles[0]), makeHaplotype(alleles[1]));
    } else {
      if (!GeneReport.isSinglePloidy(f_gene) && !isPhenotypeOnly()) {
        throw new RuntimeException("Expected two genotypes separated by a '/' but saw [" + diplotypeText + "] for " + f_gene);
      }
      diplotype = new Diplotype(f_gene, makeHaplotype(diplotypeText));
    }
    fillDiplotype(diplotype);

    return diplotype;
  }

  private boolean isPhenotypeOnly() {
    return PHENOTYPE_ONLY.contains(f_gene);
  }

  private Diplotype makeUnknownDiplotype() {
    return makeDiplotype(Diplotype.UNKNOWN);
  }

  private void fillDiplotype(Diplotype diplotype) {
    if (f_genePhenotype != null && !f_gene.startsWith("HLA")) {
      diplotype.addPhenotype(f_genePhenotype.getPhenotypeForDiplotype(diplotype));
      diplotype.addLookupKey(f_genePhenotype.getLookupKeyForDiplotype(diplotype));
    }
    if (f_gene.startsWith("HLA")) {
      diplotype.setPhenotypes(makeHlaPhenotype(diplotype));
      diplotype.setLookupKeys(makeHlaPhenotype(diplotype));
    }
  }

  /**
   * Make or retrieve a cached Haplotype object that corresponds to the given allele name
   * @param name an allele name (e.g. *20)
   * @return a Haplotype object for that allele (new or cached)
   */
  private Haplotype makeHaplotype(String name) {

    // return cache value if possible
    if (m_haplotypeCache.containsKey(name)) {
      return m_haplotypeCache.get(name);
    }

    Haplotype haplotype = new Haplotype(f_gene, name);
    if (f_genePhenotype != null) {
      haplotype.setFunction(f_genePhenotype.findHaplotypeFunction(name).orElse(UNASSIGNED_FUNCTION));
    }
    haplotype.setReference(name.equals(f_referenceAlleleName));

    m_haplotypeCache.put(name, haplotype);
    return haplotype;
  }

  private Haplotype makeHaplotype(HaplotypeMatch haplotypeMatch) {
    return makeHaplotype(haplotypeMatch.getName());
  }

  private List<String> makeHlaPhenotype(Diplotype diplotype) {
    if (!f_gene.equals("HLA-A") && !f_gene.equals("HLA-B")) {
      throw new RuntimeException("Gene not supported for HLA phenotype calling: " + f_gene);
    }
    if (diplotype.isUnknown()) {
      return new ArrayList<>();
    }
    if (f_gene.equals("HLA-A")) {
      return ImmutableList.of(diplotype.containsAllele("*31:01"));
    } else {
      List<String> phenotypes = new ArrayList<>();
      phenotypes.add(diplotype.containsAllele("*15:02"));
      phenotypes.add(diplotype.containsAllele("*57:01"));
      phenotypes.add(diplotype.containsAllele("*58:01"));
      return phenotypes;
    }
  }
}
