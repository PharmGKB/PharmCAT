package org.pharmgkb.pharmcat.reporter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSortedSet;
import org.pharmgkb.pharmcat.ParseException;
import org.pharmgkb.pharmcat.UnexpectedStateException;
import org.pharmgkb.pharmcat.definition.IncidentalFinder;
import org.pharmgkb.pharmcat.definition.model.GenePhenotype;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.model.AstrolabeCall;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;


/**
 * Factory class for spitting out {@link Diplotype} objects from known information about the gene
 *
 * This is designed to take information from either the Matcher or from Astrolabe and produce consistent results
 *
 * @author Ryan Whaley
 */
public class DiplotypeFactory {

  private final String f_gene;
  private final String f_referenceAlleleName;
  private final GenePhenotype f_genePhenotype;
  private final IncidentalFinder f_incidentalFinder;
  private final Set<Variant> f_variants;

  private Map<String,Haplotype> m_haplotypeCache = new HashMap<>();

  /**
   * public constructor
   *
   * Initialize the factory with all the necessary information to make haplotype and diplotype calls
   *
   * @param gene the gene symbol for the diplotypes this will call
   * @param variants variant calls, optional
   * @param genePhenotype a {@link GenePhenotype} object that maps haplotypes and functions to phenotypes
   */
  public DiplotypeFactory(String gene, @Nullable Collection<Variant> variants, GenePhenotype genePhenotype, IncidentalFinder incidentalFinder, @Nonnull String referenceAlleleName) {
    f_gene = gene;
    f_incidentalFinder = incidentalFinder;
    f_genePhenotype = genePhenotype;
    f_referenceAlleleName = referenceAlleleName;
    if (variants != null) {
      f_variants = ImmutableSortedSet.copyOf(variants);
    } else {
      //noinspection unchecked
      f_variants = Collections.EMPTY_SET;
    }
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

    // some genes get a second shot at the call in cases where certain variant data is known
    else {
      return ImmutableList.of();
    }
  }

  /**
   * Make diplotype objects based on AstrolabeCall objects that are parsed from astrolabe output
   */
  public List<Diplotype> makeDiplotypes(@Nonnull AstrolabeCall astrolabeCall) {
    Preconditions.checkNotNull(astrolabeCall);
    Preconditions.checkArgument(astrolabeCall.getGene().equals(f_gene));

    return makeDiplotypes(astrolabeCall.getDiplotypes());
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
    Preconditions.checkArgument(diplotypeText.contains("/"));

    String[] alleles = diplotypeText.split("/");

    Diplotype diplotype = new Diplotype(f_gene, makeHaplotype(alleles[0]), makeHaplotype(alleles[1]));
    fillDiplotype(diplotype);

    return diplotype;
  }

  private void fillDiplotype(Diplotype diplotype) {
    if (f_genePhenotype != null) {
      diplotype.setPhenotype(f_genePhenotype.makePhenotype(diplotype.printBareLookupKey()));
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
    if (f_genePhenotype != null && f_genePhenotype.getHaplotypes() != null) {
      haplotype.setGuidelineFunction(f_genePhenotype.getHaplotypes().get(name));
    }
    haplotype.setIncidental(f_incidentalFinder);
    haplotype.setReference(name.equals(f_referenceAlleleName));

    m_haplotypeCache.put(name, haplotype);
    return haplotype;
  }

  private Haplotype makeHaplotype(HaplotypeMatch haplotypeMatch) {
    String name = haplotypeMatch.getName();

    Haplotype haplotype = makeHaplotype(name);
    haplotype.setCalledFunction(haplotypeMatch.getFunction());

    return haplotype;
  }

  /**
   * Calls diplotypes based on variant data as a replacement for typical diplotype calling
   * @return a List of called Diplotypes, empty if no call
   */
  public List<Diplotype> makeOverrideDiplotypes() {
    List<Diplotype> diplotypes = new ArrayList<>();

    switch (f_gene) {
      case "SLCO1B1":
        Diplotype knownRs4149056 = callSlco1b1();
        if (knownRs4149056 != null) {
          diplotypes.add(knownRs4149056);
        }
        break;
    }

    return diplotypes;
  }

  /**
   * Calls the SLCO1B1 gene in the case when rs4149056 is available
   * @return the called Diplotype, can be null
   */
  @Nullable
  private Diplotype callSlco1b1() {
    Variant variant = f_variants.stream().filter(v -> v.getRsid() != null && v.getRsid().equals("rs4149056")).reduce((a,b) -> {
      throw new UnexpectedStateException("more than one variant found");
    }).orElse(null);

    if (variant != null && variant.getVcfCall() != null) {
      String[] alleles = variant.getVcfCall().split("[|/]");

      Diplotype dip;
      if (Arrays.equals(alleles, new String[]{"T","T"})) {
        dip = makeDiplotype("*1A/*1A");
      }
      else if (Arrays.equals(alleles, new String[]{"C","C"})) {
        dip = makeDiplotype("*5/*5");
      }
      else if ((Arrays.equals(alleles, new String[]{"T","C"})) || (Arrays.equals(alleles, new String[]{"C","T"}))) {
        dip = makeDiplotype("*1A/*5");
      }
      else {
        throw new ParseException("Unexpected genotype for " + variant);
      }

      dip.setVariant(variant);
      return dip;
    } else {
      return null;
    }
  }
}
