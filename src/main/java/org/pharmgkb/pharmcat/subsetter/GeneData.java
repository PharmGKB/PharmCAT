package org.pharmgkb.pharmcat.subsetter;

import java.util.Arrays;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;


/**
 * @author Mark Woon
 */
public class GeneData implements Comparable<GeneData> {
  public final DefinitionFile definitionFile;
  public final String gene;
  public final SortedSet<NamedAllele> haplotypes = new TreeSet<>();
  public final SortedSet<NamedAllele> extraHaplotypes = new TreeSet<>();
  public final SortedSet<NamedAllele> subsetHaplotypes = new TreeSet<>();
  public final SortedSet<NamedAllele> overlapHaplotypes = new TreeSet<>();
  public final SortedSet<NamedAllele> unusedHaplotypes = new TreeSet<>();

  public final SortedSet<VariantLocus> variants = new TreeSet<>();
  public final SortedSet<VariantLocus> missedVariants = new TreeSet<>();
  public final SortedSet<VariantLocus> unusedVariants = new TreeSet<>();


  GeneData(DefinitionFile definitionFile) {
    this.definitionFile = definitionFile;
    gene = definitionFile.getGeneSymbol();
  }

  public void calculateUnused() {
    Arrays.stream(definitionFile.getVariants())
        .filter(vl -> !variants.contains(vl) && !missedVariants.contains(vl))
        .forEach(unusedVariants::add);
    definitionFile.getNamedAlleles().stream()
        .filter(na -> !haplotypes.contains(na) && !subsetHaplotypes.contains(na) && !overlapHaplotypes.contains(na))
        .forEach(unusedHaplotypes::add);
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(gene);
  }

  @Override
  public boolean equals(Object o) {
    if (o == this) {
      return true;
    }
    if (!(o instanceof GeneData gd)) {
      return false;
    }
    return Objects.equals(gene, gd.gene);
  }

  @Override
  public int compareTo(GeneData o) {
    return gene.compareTo(o.gene);
  }
}
