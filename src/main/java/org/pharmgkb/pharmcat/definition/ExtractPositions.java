package org.pharmgkb.pharmcat.definition;

import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.collect.ImmutableSet;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.util.DataManager;
import org.pharmgkb.pharmcat.util.VcfHelper;


/**
 * This extracts the positions from the definition files and reformats them as a vcf. On the command line it takes one
 * argument:
 *
 * <ul>
 *   <li>-o output_vcf_path = A path to write the VCF file to</li>
 * </ul>
 *
 * @author Lester Carter
 * @author Ryan Whaley
 */
public class ExtractPositions {
  private static final Set<String> sf_excludedGenes = ImmutableSet.of("CYP2D6");


  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("o", "output-file", "output vcf file", true, "o")
        .addOption("a", "all-genes", "include all genes");

    cliHelper.execute(args, cli -> {
      try {
        // load the allele definition files
        DefinitionReader definitionReader = new DefinitionReader();
        definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);

        SortedSet<String> m_genes = new TreeSet<>(definitionReader.getGenes());
        if (!cli.hasOption("a")) {
          m_genes.removeAll(sf_excludedGenes);
        }
        if (m_genes.size() == 0) {
          System.err.println("Did not find any allele definitions at " + DataManager.DEFAULT_DEFINITION_DIR);
          return 1;
        }

        Path outputVcf = cli.getValidFile("o", false);
        System.out.println("Writing to " + outputVcf + "...");
        VcfHelper.extractPositions(m_genes, definitionReader, outputVcf);
        System.out.println("Done.");
        return 0;

      } catch (Exception ex) {
        ex.printStackTrace();
        return 1;
      }
    });
  }
}
