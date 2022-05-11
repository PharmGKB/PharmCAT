package org.pharmgkb.pharmcat.reporter.model.pgkb;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * PharmGKB Gene Model. This is an expanded model of the {@link AccessionObject} model that has more gene-specific
 * properties.
 */
public class Gene extends AccessionObject {
  @SerializedName("chr")
  @Expose
  private AccessionObject m_chromosome;
  @SerializedName("strand")
  @Expose
  private String m_strand;


  public AccessionObject getChromosome() {
    return m_chromosome;
  }

  public void setChromosome(AccessionObject chromosome) {
    m_chromosome = chromosome;
  }

  public String getStrand() {
    return m_strand;
  }

  public void setStrand(String strand) {
    m_strand = strand;
  }
}
