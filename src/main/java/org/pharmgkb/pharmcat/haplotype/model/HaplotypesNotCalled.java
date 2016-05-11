
package org.pharmgkb.pharmcat.haplotype.model;

import java.util.ArrayList;
import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * Haplotypes that cannot be called because of missing data in VCF file.
 */
public class HaplotypesNotCalled {

  @SerializedName("haplotype")
  @Expose
  private String haplotype;
  @SerializedName("positionsMissing")
  @Expose
  private List<PositionsMissing> positionsMissing = new ArrayList<>();

  /**
   * @return The haplotype
   */
  public String getHaplotype() {
    return haplotype;
  }

  /**
   * @param haplotype The haplotype
   */
  public void setHaplotype(String haplotype) {
    this.haplotype = haplotype;
  }

  /**
   * @return The positionsMissing
   */
  public List<PositionsMissing> getPositionsMissing() {
    return positionsMissing;
  }

  /**
   * @param positionsMissing The positionsMissing
   */
  public void setPositionsMissing(List<PositionsMissing> positionsMissing) {
    this.positionsMissing = positionsMissing;
  }

}
