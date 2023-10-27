package org.pharmgkb.pharmcat.haplotype;

import java.util.Collection;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Metadata;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;


/**
 * This class collects results for the {@link NamedAlleleMatcher}.
 *
 * @author Mark Woon
 */
public class ResultBuilder {
  private final DefinitionReader m_definitionReader;
  private final Result m_result = new Result();
  private final boolean m_topCandidatesOnly;
  private final boolean m_findCombinations;
  private final boolean m_callCyp2d6;


  public ResultBuilder(DefinitionReader definitionReader, boolean topCandidatesOnly, boolean findCombinations,
      boolean callCyp2d6) {
    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
    m_topCandidatesOnly = topCandidatesOnly;
    m_findCombinations = findCombinations;
    m_callCyp2d6 = callCyp2d6;
  }

  public Result build() {
    return m_result;
  }


  public ResultBuilder forFile(VcfFile vcfFile, Map<String, Collection<String>> warnings) {
    Preconditions.checkNotNull(vcfFile);

    m_result.setMetadata(new Metadata(NamedAlleleMatcher.VERSION, m_definitionReader.getGenomeBuild(),
        PathUtils.getFilename(vcfFile.getFile()), new Date(), m_topCandidatesOnly, m_findCombinations, m_callCyp2d6));
    if (warnings != null) {
      m_result.setVcfWarnings(warnings);
    }
    return this;
  }


  /**
   * Builds the result for gene when VCF has with no samples for it.
   */
  protected ResultBuilder gene(String gene, MatchData matchData) {
    Preconditions.checkNotNull(gene);
    m_result.addGeneCall(initGeneCall(gene, matchData, null));
    return this;
  }


  /**
   * Adds diplotype results for specified gene.
   */
  protected ResultBuilder diplotypes(String gene, MatchData matchData, List<DiplotypeMatch> matches) {
    Preconditions.checkNotNull(gene);
    return diplotypes(gene, matchData, matches, null);
  }

  /**
   * Adds diplotype results for specified gene.
   */
  protected ResultBuilder diplotypes(String gene, MatchData matchData, List<DiplotypeMatch> matches,
      @Nullable List<MessageAnnotation> warnings) {
    Preconditions.checkNotNull(gene);

    GeneCall geneCall = initGeneCall(gene, matchData, warnings);
    // get haplotype/diplotype info
    for (DiplotypeMatch dm : matches) {
      geneCall.addDiplotype(dm);
    }

    m_result.addGeneCall(geneCall);
    return this;
  }


  /**
   * Add haplotype results for specified gene.
   * <p>
   * This should only be used when we can't get diplotypes but still need to track potential haplotypes (e.g. DPYD).
   */
  protected ResultBuilder haplotypes(String gene, MatchData matchData, List<HaplotypeMatch> matches) {
    return haplotypes(gene, matchData, matches, null);
  }

  /**
   * Add haplotype results for specified gene.
   * <p>
   * This should only be used when we can't get diplotypes but still need to track potential haplotypes (e.g. DPYD).
   */
  protected ResultBuilder haplotypes(String gene, MatchData matchData, List<HaplotypeMatch> matches,
      @Nullable List<MessageAnnotation> warnings) {
    Preconditions.checkNotNull(gene);

    GeneCall geneCall = initGeneCall(gene, matchData, warnings);
    if (matches != null) {
      geneCall.addHaplotypeMatches(matches);
    }

    m_result.addGeneCall(geneCall);
    return this;
  }


  private GeneCall initGeneCall(String gene, MatchData matchData, @Nullable List<MessageAnnotation> warnings) {
    Set<String> matchableHaps = matchData.getHaplotypes().stream()
        .map(NamedAllele::getName)
        .collect(Collectors.toSet());
    Set<String> uncallableHaplotypes = m_definitionReader.getHaplotypes(gene).stream()
        .map(NamedAllele::getName)
        .filter(n -> !matchableHaps.contains(n))
        .collect(Collectors.toSet());
    Set<String> ignoredHaplotypes;
    DefinitionExemption exemption = m_definitionReader.getExemption(gene);
    if (exemption != null) {
      uncallableHaplotypes = uncallableHaplotypes.stream()
          .filter(h -> !exemption.shouldIgnoreAllele(h))
          .collect(Collectors.toSet());
      ignoredHaplotypes = exemption.getIgnoredAlleles().stream()
          .map(String::toUpperCase)
          .collect(Collectors.toSet());
    } else {
      ignoredHaplotypes = new HashSet<>();
    }

    DefinitionFile definitionFile = m_definitionReader.getDefinitionFile(gene);
    GeneCall geneCall = new GeneCall(definitionFile.getSource(), definitionFile.getVersion(),
        definitionFile.getChromosome(), gene, matchData, uncallableHaplotypes, ignoredHaplotypes, warnings);

    // get position info
    for (VariantLocus variant : matchData.getPositions()) {
      geneCall.add(new Variant(variant, matchData.getSampleAllele(variant.getPosition())));
    }
    geneCall.finalizeVariants();

    return geneCall;
  }
}
