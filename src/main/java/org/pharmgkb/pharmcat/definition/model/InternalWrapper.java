package org.pharmgkb.pharmcat.definition.model;

import java.io.IOException;
import java.util.Arrays;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import org.pharmgkb.pharmcat.util.VcfHelper;


/**
 * This wrapper hides methods that should only be used by PharmCAT internally.
 *
 * @author Mark Woon
 */
public class InternalWrapper {

  /**
   * Private constructor for utility class.
   */
  private InternalWrapper() {
  }


  /**
   * Filters out structural variant alleles, they have no definition to match against.
   * <p>
   * This should only be used during data ingestion/prep.
   */
  public static void removeStructuralVariants(DefinitionFile definitionFile) {
    definitionFile.removeStructuralVariants();
  }


  /**
   * Removes ignored allele specified in {@link DefinitionExemption}.
   * <p>
   * This should only be used during data ingestion/prep.
   */
  public static void removeIgnoredNamedAlleles(DefinitionFile definitionFile, DefinitionExemption exemption) {
    System.out.println("  Removing " + exemption.getIgnoredAlleles());
    definitionFile.removeIgnoredNamedAlleles(exemption);
  }


  /**
   * Removes any unused positions and ignored positions specified in {@link DefinitionExemption}.
   * <p>
   * This should only be used during data ingestion/prep.
   */
  public static void removeIgnoredPositions(DefinitionFile definitionFile, DefinitionExemption exemption) {
    SortedSet<VariantLocus> ignoredPositions = Arrays.stream(definitionFile.getVariants())
        .filter(exemption::shouldIgnorePosition)
        .collect(Collectors.toCollection(TreeSet::new));
    if (exemption.getIgnoredPositions().size() != ignoredPositions.size()) {
      throw new IllegalStateException("Should have " + exemption.getIgnoredPositions().size() + " ignored positions, " +
          "but only found " + ignoredPositions.size());
    }
    definitionFile.removeIgnoredPositions(ignoredPositions, true);
  }

  /**
   * Removes any unused positions and the specified {@code ignoredPositions}.
   * <p>
   * This should only be used during data ingestion/prep.
   */
  public static void removeIgnoredPositions(DefinitionFile definitionFile, SortedSet<VariantLocus> ignoredPositions,
      boolean verbose) {
    definitionFile.removeIgnoredPositions(ignoredPositions, verbose);
  }


  /**
   * Translate variants from CPIC to VCF (i.e. {@code cpicAlleles} to {@code alleles}).
   * <p>
   * This should only be used during data ingestion/prep.
   */
  public static void doVcfTranslation(DefinitionFile definitionFile, VcfHelper vcfHelper) throws IOException {
    definitionFile.doVcfTranslation(vcfHelper);
  }


  /**
   * Resets this {@link DefinitionFile}'s named alleles to the specified set.
   * <p>
   * This should only be used during data ingestion/prep.
   */
  public static void resetNamedAlleles(DefinitionFile definitionFile, SortedSet<NamedAllele> namedAlleles) {
    definitionFile.resetNamedAlleles(namedAlleles);
  }

  /**
   * Adds a {@link NamedAllele} to this {@link DefinitionFile}.
   * <p>
   * This should only be used during data ingestion/prep.
   */
  public static void addNamedAllele(DefinitionFile definitionFile, NamedAllele namedAllele) {
    definitionFile.addNamedAllele(namedAllele);
  }
}
