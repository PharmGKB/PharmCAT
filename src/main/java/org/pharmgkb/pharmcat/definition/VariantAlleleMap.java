package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Comparator;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import com.google.gson.Gson;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * Class to map genomic loci to the alleles they help define
 *
 * @author Ryan Whaley
 */
public class VariantAlleleMap {

  private Multimap<Integer,String> m_variantAlleleMap = TreeMultimap.create(
      Comparator.naturalOrder(),
      HaplotypeNameComparator.getComparator()
  );

  /**
   * public constructor. will create the map for the given gene (by symbol)
   * @param gene the gene symbol (all caps)
   */
  public VariantAlleleMap(@Nonnull String gene) throws IOException {
    Gson gson = new Gson();

    Path definitionPath = DataManager.DEFAULT_DEFINITION_DIR.resolve(gene+"_translation.json");

    if (!Files.isRegularFile(definitionPath)) {
      return;
    }

    DefinitionFile definitionFile;
    try (BufferedReader reader = Files.newBufferedReader(definitionPath)) {
      definitionFile = gson.fromJson(reader, DefinitionFile.class);
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    if (definitionFile == null) {
      return;
    }

    VariantLocus[] allVariants = definitionFile.getVariants();

    // purposely skip the first allele since it's "default"
    for (int j = 1; j < definitionFile.getNamedAlleles().size(); j++) {
      NamedAllele namedAllele = definitionFile.getNamedAlleles().get(j);

      String[] definingAlleles = namedAllele.getAlleles();
      for (int i = 0; i < definingAlleles.length; i++) {
        String alleleValue = definingAlleles[i];
        VariantLocus locus = allVariants[i];

        if (alleleValue != null) {
          m_variantAlleleMap.put(locus.getPosition(), namedAllele.getName());
        }
      }
    }
  }

  /**
   * Gets the allele names (e.g. *2) that the given chromosomal position for the gene is used to define
   * @param position a chromosomal position
   * @return a Collection of allele names
   */
  @Nullable
  public Collection<String> getAlleles(int position) {
    return m_variantAlleleMap.get(position);
  }
}
