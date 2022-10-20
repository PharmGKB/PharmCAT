package org.pharmgkb.pharmcat.definition;

import java.util.HashMap;
import java.util.Map;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.StringUtils;


/**
 * A map that stores the name of the reference allele for each gene.
 */
public class ReferenceAlleleMap {
  private final Map<String,String> f_refAlleleForGene = new HashMap<>();


  ReferenceAlleleMap(DefinitionReader definitionReader) {
    f_refAlleleForGene.put("MT-RNR1", "Reference");
    // use reference prop
    definitionReader.getGenes()
        .forEach(g -> f_refAlleleForGene.put(g, definitionReader.getHaplotypes(g).first().getName()));
  }


  public String get(String geneSymbol) {
    Preconditions.checkArgument(StringUtils.isNotBlank(geneSymbol), "Must supply a gene symbol");
    return f_refAlleleForGene.get(geneSymbol);
  }
}
