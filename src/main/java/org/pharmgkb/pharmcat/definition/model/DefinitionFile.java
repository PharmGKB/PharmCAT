package org.pharmgkb.pharmcat.definition.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.gson.annotations.SerializedName;


/**
 * This class represents one complete allele translation set for a gene.
 *
 * @author Ryan Whaley
 */
public class DefinitionFile {
  public static final String FORMAT_VERSION = "1";
  // metadata
  @SerializedName("formatVersion")
  private String m_formatVersion = FORMAT_VERSION;
  @SerializedName("contentVersion")
  private String m_contentVersion;
  @SerializedName("modificationDate")
  private Date m_modificationDate;
  @SerializedName("gene")
  private String m_geneSymbol;
  @SerializedName("orientation")
  private String m_orientation;
  @SerializedName("chromosome")
  private String m_chromosome;
  @SerializedName("genomeBuild")
  private String m_genomeBuild;
  @SerializedName("refSeqChromosomeId")
  private String m_refSeqChromosome;
  @SerializedName("refSeqGeneId")
  private String m_refSeqGene;
  @SerializedName("refSeqProteinId")
  private String m_refSeqProtein;
  @SerializedName("notes")
  private List<String> m_notes;
  @SerializedName("populations")
  private SortedSet<String> m_populations;
  // definitions
  @SerializedName("variants")
  private VariantLocus[] m_variants;
  @SerializedName("namedAlleles")
  private List<NamedAllele> m_namedAlleles;


  /**
   * The format version of the definition file.
   */
  public String getFormatVersion() {
    return m_formatVersion;
  }

  public void setFormatVersion(String formatVersion) {
    m_formatVersion = formatVersion;
  }

  /**
   * The version of the actual data content, this will iterate when curators have changed information in any of the
   * member named alleles.
   */
  public String getContentVersion() {
    return m_contentVersion;
  }

  public void setContentVersion(String contentVersion) {
    m_contentVersion = contentVersion;
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
   * The symbol for the gene these alleles are on
   */
  public String getGeneSymbol() {
    return m_geneSymbol;
  }

  public void setGeneSymbol(String geneSymbol) {
    m_geneSymbol = geneSymbol;
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
   * The name of the chromosome this translation is on.
   */
  public String getChromosome() {
    return m_chromosome;
  }

  public void setChromosome(String chromosome) {
    m_chromosome = chromosome;
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
   * The RefSeq identifier for the chromosome this translation is on (should agree with build).
   */
  public String getRefSeqChromosome() {
    return m_refSeqChromosome;
  }

  public void setRefSeqChromosome(String refSeqChromosome) {
    m_refSeqChromosome = refSeqChromosome;
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
   * The RefSeq identifier for the protein the gene translates to
   */
  public String getRefSeqProtein() {
    return m_refSeqProtein;
  }

  public void setRefSeqProtein(String refSeqProtein) {
    m_refSeqProtein = refSeqProtein;
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
   * The different populations used for variant frequencies in this translation.
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
   * The {@link VariantLocus} objects used to define alleles in this translation
   */
  public VariantLocus[] getVariants() {
    return m_variants;
  }

  public void setVariants(VariantLocus[] variants) {
    m_variants = variants;
  }


  /**
   * All the named alleles defined in this translation
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


  @Override
  public String toString() {
    return "Allele definition for " + m_geneSymbol + " (v." + m_contentVersion + ")";
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (!(o instanceof DefinitionFile)) return false;
    DefinitionFile that = (DefinitionFile)o;
    return Objects.equals(m_formatVersion, that.getFormatVersion()) &&
        Objects.equals(m_contentVersion, that.getContentVersion()) &&
        Objects.equals(m_modificationDate, that.getModificationDate()) &&
        Objects.equals(m_geneSymbol, that.getGeneSymbol()) &&
        Objects.equals(m_orientation, that.getOrientation()) &&
        Objects.equals(m_chromosome, that.getChromosome()) &&
        Objects.equals(m_genomeBuild, that.getGenomeBuild()) &&
        Objects.equals(m_refSeqChromosome, that.getRefSeqChromosome()) &&
        Objects.equals(m_refSeqGene, that.getRefSeqGene()) &&
        Objects.equals(m_refSeqProtein, that.getRefSeqProtein()) &&
        Objects.equals(m_notes, that.getNotes()) &&
        Objects.equals(m_populations, that.getPopulations()) &&
        Arrays.equals(m_variants, that.getVariants()) &&
        Objects.equals(m_namedAlleles, that.getNamedAlleles());
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_formatVersion, m_contentVersion, m_modificationDate, m_geneSymbol, m_orientation,
        m_chromosome, m_genomeBuild, m_refSeqChromosome, m_refSeqGene, m_refSeqProtein, m_notes, m_populations,
        m_variants, m_namedAlleles);
  }
}
