package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * A map that stores the name of the reference allele for each gene
 */
public class ReferenceAlleleMap {
  private final Map<String,String> f_refAlleleForGene = new HashMap<>();

  public ReferenceAlleleMap() {
    f_refAlleleForGene.put("CYP2D6", "*1");
    f_refAlleleForGene.put("G6PD", "B (wildtype)");
    f_refAlleleForGene.put("MT-RNR1", "Reference");
    try {
      DefinitionReader definitionReader = new DefinitionReader();
      definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
      definitionReader.getGenes()
          .forEach(g -> f_refAlleleForGene.put(g, definitionReader.getHaplotypes(g).first().getName()));
    } catch (IOException ex) {
      throw new RuntimeException("Error loading reference allele map", ex);
    }
  }

  public String get(String geneSymbol) {
    Preconditions.checkArgument(StringUtils.isNotBlank(geneSymbol), "Must supply a gene symbol");

    if (f_refAlleleForGene.containsKey(geneSymbol)) {
      return f_refAlleleForGene.get(geneSymbol);
    } else {
      throw new RuntimeException("No reference allele for gene [" + geneSymbol + "]");
    }
  }
}
