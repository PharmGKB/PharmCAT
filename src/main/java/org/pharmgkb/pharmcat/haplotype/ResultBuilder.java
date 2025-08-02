package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.UnphasedDiplotypePriority;
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
  private Path m_sampleMetadataFile;


  public ResultBuilder(DefinitionReader definitionReader, boolean topCandidatesOnly, boolean findCombinations,
      boolean callCyp2d6) {
    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
    m_topCandidatesOnly = topCandidatesOnly;
    m_findCombinations = findCombinations;
    m_callCyp2d6 = callCyp2d6;
  }

  public Result build(Env env) throws IOException {
    if (m_sampleMetadataFile != null) {
      Metadata metadata = m_result.getMetadata();
      Map<String, String> sampleData = env.getSampleMetadata(m_sampleMetadataFile, metadata.getSampleId(), true);
      if (sampleData != null && !sampleData.isEmpty()) {
        metadata.setSampleProps(sampleData);
      }
    }
    for (GeneCall call : m_result.getGeneCalls()) {
      if (call.getDiplotypes().size() > 1 && !call.isEffectivelyPhased()) {
        DefinitionExemption exemption = m_definitionReader.getExemption(call.getGene());
        if (exemption != null && exemption.hasUnphasedDiplotypePriorities()) {
          Set<String> dips = call.getDiplotypes().stream()
              .map(DiplotypeMatch::getName)
              .collect(Collectors.toCollection(TreeSet::new));
          for (UnphasedDiplotypePriority udp : exemption.getUnphasedDiplotypePriorities()) {
            if (dips.containsAll(udp.getList())) {
              DiplotypeMatch priorityDip = call.getDiplotypes().stream()
                  .filter(dm -> dm.getName().equals(udp.getPick()))
                  .findAny()
                  // this should never happen
                  .orElseThrow(() -> new IllegalStateException("Cannot find " + udp.getPick()));
              call.setPriorityDiplotype(priorityDip);
              call.addWarning(new MessageAnnotation(MessageAnnotation.TYPE_NOTE,
                  "unphased-priority",
                  "Unphased " + call.getGene() + " variants resulted in multiple calls.  " +
                      "PharmCAT is picking a single call based on frequency data.  " +
                      "Please consult the documentation for details."));
              break;
            }
          }
        }
      }
    }
    return m_result;
  }


  public ResultBuilder forFile(VcfFile vcfFile, Map<String, Collection<String>> warnings, String sampleId,
      @Nullable Path sampleMetadataFile) {
    Preconditions.checkNotNull(vcfFile);

    Metadata metadata = new Metadata(NamedAlleleMatcher.VERSION, m_definitionReader.getGenomeBuild(),
        PathUtils.getFilename(vcfFile.getFile()), new Date(), m_topCandidatesOnly, m_findCombinations, m_callCyp2d6,
        sampleId);
    m_result.setMetadata(metadata);
    if (warnings != null) {
      m_result.setVcfWarnings(warnings);
    }
    m_sampleMetadataFile = sampleMetadataFile;
    return this;
  }


  /**
   * Builds the result for gene when VCF has no samples for it.
   */
  protected ResultBuilder noCall(String gene, MatchData matchData) {
    Preconditions.checkNotNull(gene);
    m_result.addGeneCall(initGeneCall(gene, matchData, null));
    return this;
  }

  /**
   * Builds the result for gene when VCF has no call can be made.
   */
  protected ResultBuilder noCall(String gene, MatchData matchData, @Nullable List<MessageAnnotation> warnings) {
    Preconditions.checkNotNull(gene);
    m_result.addGeneCall(initGeneCall(gene, matchData, warnings));
    return this;
  }


  /**
   * Adds diplotype results for a specified gene.
   */
  protected ResultBuilder diplotypes(String gene, MatchData matchData, SortedSet<DiplotypeMatch> matches) {
    Preconditions.checkNotNull(gene);
    return diplotypes(gene, matchData, matches, null);
  }

  /**
   * Adds diplotype results for a specified gene.
   */
  protected ResultBuilder diplotypes(String gene, MatchData matchData, SortedSet<DiplotypeMatch> matches,
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
   * Add haplotype results for a specified gene.
   * <p>
   * This should only be used when we can't get diplotypes but still need to track potential haplotypes (e.g. DPYD).
   */
  protected ResultBuilder haplotypes(String gene, MatchData matchData, List<HaplotypeMatch> matches) {
    return haplotypes(gene, matchData, matches, null);
  }

  /**
   * Add haplotype results for a specified gene.
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

    if (!matchData.getMissingAmp1Positions().isEmpty()) {
      StringBuilder builder = new StringBuilder("Missing variant");
      if (matchData.getMissingAmp1Positions().size() > 1) {
        builder.append("s");
      }
      builder.append(" required to meet AMP Tier 1 requirements:  ")
          .append(String.join(", ", matchData.getMissingRequiredPositions()))
          .append(". See https://www.clinpgx.org/ampAllelesToTest for details.");
      if (warnings == null) {
        warnings = new ArrayList<>();
      }
      warnings.add(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "missing-amp1-position", builder.toString()));
    }

    DefinitionFile definitionFile = m_definitionReader.getDefinitionFile(gene);
    GeneCall geneCall = new GeneCall(definitionFile.getSource(), definitionFile.getVersion(),
        definitionFile.getChromosome(), gene, matchData, uncallableHaplotypes, warnings);

    // get position info
    for (VariantLocus variant : matchData.getPositions()) {
      geneCall.add(new Variant(variant, matchData.getSampleAllele(variant.getPosition())));
    }
    geneCall.finalizeVariants();

    return geneCall;
  }
}
