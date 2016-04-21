package org.pharmgkb.pharmcat.definition.model;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This class represents one complete allele translation set for a gene.
 *
 * @author Ryan Whaley
 */
public class AlleleTranslation {

  private String m_formatVersion;
  private String m_contentVersion;
  private String m_geneSymbol;
  private String m_refSeqGene;
  private String m_orientation;
  private Date m_modificationDate;
  private String m_genomeBuild;
  private String m_chromoName;
  private String m_refSeqChromo;
  private String m_refSeqProtein;
  private VariantLocus[] m_variants;
  private List<String> m_notes;
  private SortedSet<String> m_populations;
  private List<NamedAllele> m_namedAlleles;


  /**
   * The format version of the original allele sheet
   */
  public String getFormatVersion() {
    return m_formatVersion;
  }

  public void setFormatVersion(String formatVersion) {
    m_formatVersion = formatVersion;
  }

  /**
   * The version of the actual data content, this will iterate when curators have changed information in any of the
   * member named alleles
   */
  public String getContentVersion() {
    return m_contentVersion;
  }

  public void setContentVersion(String contentVersion) {
    m_contentVersion = contentVersion;
  }

  /**
   * The symbol for the gene these alleles are on
   */
  public String getGeneSymbol() {
    return m_geneSymbol;
  }

  public void setGeneSymbol(String geneSymbol) {
    m_geneSymbol = geneSymbol;
  }

  /**
   * The RefSeq identifier for the gene in this translation
   */
  public String getRefSeqGene() {
    return m_refSeqGene;
  }

  public void setRefSeqGene(String refSeqGene) {
    m_refSeqGene = refSeqGene;
  }

  /**
   * The orientation of the gene relative to the chromosome
   */
  public String getOrientation() {
    return m_orientation;
  }

  public void setOrientation(String orientation) {
    m_orientation = orientation;
  }

  /**
   * The date this file should be considered to be last modified (should be manually set by curators)
   */
  public Date getModificationDate() {
    return m_modificationDate;
  }

  public void setModificationDate(Date modificationDate) {
    m_modificationDate = modificationDate;
  }

  /**
   * The human genome assembly (build) the positions in this translation are from (e.g. b38 or b37)
   */
  public String getGenomeBuild() {
    return m_genomeBuild;
  }

  public void setGenomeBuild(String genomeBuild) {
    m_genomeBuild = genomeBuild;
  }

  /**
   * The name of the chromosome this translation is on
   */
  public String getChromoName() {
    return m_chromoName;
  }

  public void setChromoName(String chromoName) {
    m_chromoName = chromoName;
  }

  /**
   * The RefSeq identifier for the chromosome this translation is on (should agree with build)
   */
  public String getRefSeqChromo() {
    return m_refSeqChromo;
  }

  public void setRefSeqChromo(String refSeqChromo) {
    m_refSeqChromo = refSeqChromo;
  }

  /**
   * The RefSeq identifier for the protein the gene translates to
   */
  public String getRefSeqProtein() {
    return m_refSeqProtein;
  }

  public void setRefSeqProtein(String refSeqProtein) {
    m_refSeqProtein = refSeqProtein;
  }

  /**
   * The {@link VariantLocus} objects used to define alleles in this translation
   */
  public VariantLocus[] getVariants() {
    return m_variants;
  }

  public void setVariants(VariantLocus[] variants) {
    m_variants = variants;
  }

  /**
   * The general notes that are relevant to this translation
   */
  public List<String> getNotes() {
    return m_notes;
  }

  public void setNotes(List<String> notes) {
    m_notes = notes;
  }

  public void addNote(String note) {
    if (m_notes == null) {
      m_notes = new ArrayList<>();
    }
    m_notes.add(note);
  }

  /**
   * The different populations used for variant frequencies in this translation
   */
  public SortedSet<String> getPopulations() {
    return m_populations;
  }

  public void setPopulations(SortedSet<String> populations) {
    m_populations = populations;
  }

  public void addPopulation(String population) {
    if (m_populations == null) {
      m_populations = new TreeSet<>();
    }
    m_populations.add(population);
  }

  /**
   * All the named alleles defined in this translation
   * @return
   */
  public List<NamedAllele> getNamedAlleles() {
    return m_namedAlleles;
  }

  public void setNamedAlleles(List<NamedAllele> namedAlleles) {
    m_namedAlleles = namedAlleles;
  }

  public void addNamedAllele(NamedAllele allele) {
    if (m_namedAlleles == null) {
      m_namedAlleles = new ArrayList<>();
    }
    m_namedAlleles.add(allele);
  }
}
