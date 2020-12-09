package org.pharmgkb.pharmcat.definition.model;

import java.util.*;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.Nullable;


/**
 * This class represents one complete allele translation set for a gene.
 *
 * @author Ryan Whaley
 */
public class DefinitionFile {
  public static final String FORMAT_VERSION = "1";
  // metadata
  @Expose
  @SerializedName("formatVersion")
  private String m_formatVersion = FORMAT_VERSION;
  @Expose
  @SerializedName("modificationDate")
  private Date m_modificationDate;
  @Expose
  @SerializedName("gene")
  private String m_geneSymbol;
  @Expose
  @SerializedName("orientation")
  private String m_orientation;
  @Expose
  @SerializedName("chromosome")
  private String m_chromosome;
  @Expose
  @SerializedName("genomeBuild")
  private String m_genomeBuild;
  @Expose
  @SerializedName("refSeqChromosomeId")
  private String m_refSeqChromosome;
  @Expose
  @SerializedName("refSeqGeneId")
  private String m_refSeqGene;
  @Expose
  @SerializedName("refSeqProteinId")
  private String m_refSeqProtein;
  @Expose
  @SerializedName("notes")
  private List<String> m_notes;
  @Expose
  @SerializedName("populations")
  private SortedSet<String> m_populations;
  // definitions
  @Expose
  @SerializedName("variants")
  private VariantLocus[] m_variants;
  @Expose
  @SerializedName("variantAlleles")
  private List<Set<String>> m_variantAlleles;
  private Map<VariantLocus, Set<String>> m_variantAllelesMap;
  @Expose
  @SerializedName("namedAlleles")
  private List<NamedAllele> m_namedAlleles;
  private final SortedMap<String, VariantLocus> m_rsidMap = new TreeMap<>();


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
   * The {@link VariantLocus} objects used to define {@link NamedAllele}s in this translation
   */
  public VariantLocus[] getVariants() {
    return m_variants;
  }

  public void setVariants(VariantLocus[] variants) {
    m_variants = variants;
    for (VariantLocus varLoc : variants) {
      if (varLoc.getRsid() != null) {
        m_rsidMap.put(varLoc.getRsid(), varLoc);
      }
    }
  }

  public @Nullable VariantLocus getVariantByRsid(String rsid) {
    Preconditions.checkNotNull(rsid);
    Preconditions.checkArgument(rsid.startsWith("rs"));
    return m_rsidMap.get(rsid);
  }


  /**
   * The expected alleles for {@link VariantLocus}'s used to define {@link NamedAllele}s in this translation.
   *
   * @throws IllegalStateException if m_variantAlleles have not been generated (using {@link #generateVariantAlleles()}
   */
  public List<Set<String>> getVariantAlleles() {
    return m_variantAlleles;
  }

  public Set<String> getVariantAlleles(VariantLocus vl) {
    if (m_variantAllelesMap == null) {
      if (m_variantAlleles == null) {
        throw new IllegalStateException("Variant alleles have not been generated yet.");
      }
      if (m_variantAlleles.size() < m_variants.length) {
        throw new IllegalStateException("More variants than there are alleles.");
      }
      m_variantAllelesMap = new HashMap<>();
      for (int x = 0; x < m_variants.length; x += 1) {
        m_variantAllelesMap.put(m_variants[x], m_variantAlleles.get(x));
      }
    }
    return m_variantAllelesMap.get(vl);
  }

  /**
   * Generates variant alleles.  Must be called before using {@link #getVariantAlleles()}.
   */
  public void generateVariantAlleles() {

    m_variantAlleles = new ArrayList<>();
    for (VariantLocus varLoc : m_variants) {
      m_variantAlleles.add(
          m_namedAlleles.stream()
              .map(na -> na.getAllele(varLoc))
              .filter(Objects::nonNull)
              .collect(Collectors.toSet())
      );
    }
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
    return "Allele definition for " + m_geneSymbol;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) return true;
    if (!(o instanceof DefinitionFile)) return false;
    DefinitionFile that = (DefinitionFile)o;
    return Objects.equals(m_formatVersion, that.getFormatVersion()) &&
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
        Objects.equals(m_variantAlleles, that.getVariantAlleles()) &&
        Objects.equals(m_namedAlleles, that.getNamedAlleles());
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_formatVersion, m_modificationDate, m_geneSymbol, m_orientation,
        m_chromosome, m_genomeBuild, m_refSeqChromosome, m_refSeqGene, m_refSeqProtein, m_notes, m_populations,
        m_variants, m_namedAlleles);
  }
}
