package org.pharmgkb.pharmcat.haplotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;


/**
 * Read in haplotype definition files.
 *
 * @author Mark Woon
 */
public class DefinitionReader {
  private ListMultimap<String, Variant> m_haplotypePositions = ArrayListMultimap.create();
  private ListMultimap<String, Haplotype> m_haplotypes = ArrayListMultimap.create();
  private Map<String, DefinitionFile> m_files = new HashMap<>();


  public Map<String, DefinitionFile> getDefinitionFiles() {
    return m_files;
  }

  public ListMultimap<String, Variant> getHaplotypePositions() {
    return m_haplotypePositions;
  }

  public ListMultimap<String, Haplotype> getHaplotypes() {
    return m_haplotypes;
  }


  public void read(Path path) throws IOException {

    if (Files.isDirectory(path)) {
      Files.list(path)
          .filter(f -> f.toString().endsWith(".tsv"))
          .forEach(this::readFile);
    } else {
      readFile(path);
    }
  }


  private void readFile(@Nonnull Path file) {

    Preconditions.checkNotNull(file);
    Preconditions.checkArgument(Files.isRegularFile(file));
    System.out.println(file);

    try (BufferedReader bufferedReader = Files.newBufferedReader(file)) {

      DefinitionFile inProcessFile = new DefinitionFile(file.toString());
      ArrayList<Variant> variants = new ArrayList<>();
      ArrayList<Haplotype> haplotypes = new ArrayList<>();

      String line;
      boolean NormalFunctionAllele = true;

      while ((line = bufferedReader.readLine()) != null) {
        //System.out.println(line);

        String[] fields = line.split("\t");
        for (int i = 0; i < fields.length; i++) {
          fields[i] = fields[i].trim();
        }

        if (fields[0].equals("FormatVersion")) {
          if (fields.length > 1) {
            inProcessFile.setFormatVersion(fields[1]);
          } else {
            inProcessFile.setFormatVersion("");
          }

        } else if (fields[0].equals("GeneName")) {
          if (fields.length > 1) {
            inProcessFile.setGeneName(fields[1]);
          } else {
            inProcessFile.setGeneName("");
          }
          inProcessFile.setGeneName(fields[1]);
        } else if (fields[0].equals("GeneRefSeq")) {
          inProcessFile.setGeneID(fields[1]);
        } else if (fields[0].equals("GeneOrientation")) {
          if (fields.length > 1) {
            inProcessFile.setGeneOrientation(fields[1]);
          } else {
            inProcessFile.setGeneOrientation("");
          }
        } else if (fields[0].equals("ContentDate")) {
          if (fields.length > 1) {
            inProcessFile.setContentDate(fields[1]);
          } else {
            inProcessFile.setContentDate("");
          }
        } else if (fields[0].equals("ContentVersion")) {
          if (fields.length > 1) {
            inProcessFile.setContentVersion(fields[1]);
          } else {
            inProcessFile.setContentVersion("");
          }
        } else if (fields[0].equals("GenomeBuild")) {
          if (fields.length > 1) {
            inProcessFile.setGenomeBuild(fields[1]);
          } else {
            inProcessFile.setGenomeBuild("");
          }
        } else if (fields[0].equals("ChrName")) {
          if (fields.length > 1) {
            inProcessFile.setChromosome(fields[1]);
          } else {
            inProcessFile.setChromosome("");
          }
        } else if (fields[0].equals("ChrRefSeq")) {
          if (fields.length > 1) {
            inProcessFile.setChromosomeID(fields[1]);
          } else {
            inProcessFile.setChromosomeID("");
          }
        } else if (fields[0].equals("ProteinRefSeq")) {
          if (fields.length > 1) {
            inProcessFile.setProteinID(fields[1]);
          } else {
            inProcessFile.setProteinID("");
          }


        } else if (fields[0].equals("NumVariants")) {
          inProcessFile.setNumberOfVariants(fields[1]);
          for (int i = 0; i < inProcessFile.getNumberOfVariants(); i++) {
            Variant newVariant = new Variant(inProcessFile.getChromosome(), inProcessFile.getGeneName());
            variants.add(newVariant);
          }
        } else if (fields[0].equals("ResourceNote")) {
          if (fields.length >= 4) {
            for (int i = 4; i < fields.length; i++) {
              variants.get(i - 4).addProteinEffect(fields[i]);
            }
          }

        } else if (fields[0].equals("ProteinNote")) {
          if (fields.length >= 4) {
            for (int i = 4; i < fields.length; i++) {
              variants.get(i - 4).addProteinEffect(fields[i]);
            }
          }
        } else if (fields[0].equals("ChrPosition")) {
          for (int i = 4; i < fields.length; i++) {
            variants.get(i - 4).addHGVSg(fields[i]);
            variants.get(i - 4).setStartPOS();

          }
        } else if (fields[0].equals("GenePosition")) {
          for (int i = 4; i < fields.length; i++) {
            variants.get(i - 4).setGenePosition(fields[i]);
          }
        } else if (fields[0].equals("rsID")) {
          if (fields.length >= 4) {
            for (int i = 4; i < fields.length; i++) {
              variants.get(i - 4).setRsid(fields[i]);
            }
          }

        } else if (fields[0].equals("Allele")) {
          ArrayList<String> alleles = new ArrayList<>();
          ArrayList<Variant> hapVariants = new ArrayList<>();

          int forLength = variants.size() + 4;
          if (forLength > fields.length) {
            forLength = fields.length;
          }


          for (int i = 4; i < forLength; i++) {

            if (NormalFunctionAllele) {
              variants.get(i - 4).setREF(fields[i]);
              alleles.add(fields[i]);
              hapVariants.add(variants.get(i - 4));

            } else {
              if (!fields[i].trim().equals("")) {
                variants.get(i - 4).addALT(fields[i]);
                alleles.add(fields[i]);
                hapVariants.add(variants.get(i - 4));
              }
            }

          }
          NormalFunctionAllele = false;
          haplotypes.add(new Haplotype(fields[1], fields[2], hapVariants, fields[3], alleles));

        }
      }
      for (int i = 0; i < haplotypes.size(); i++) {
        m_haplotypes.put(inProcessFile.getGeneName(), haplotypes.get(i));
      }
      for (int i = 0; i < variants.size(); i++) {
        m_haplotypePositions.put(inProcessFile.getGeneName(), variants.get(i));
      }

      m_files.put(inProcessFile.getGeneName(), inProcessFile);


    } catch (Exception ex) {
      throw new RuntimeException("Failed to parse " + file, ex);
    }
  }
}
