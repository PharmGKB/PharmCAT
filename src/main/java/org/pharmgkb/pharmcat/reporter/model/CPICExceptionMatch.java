package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.List;
import javax.annotation.Nonnull;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * @author Ryan Whaley
 */
public class CPICExceptionMatch {
  
  @Expose
  @SerializedName("genes")
  private List<String> m_genes = new ArrayList<>();
  @Expose
  @SerializedName("hapsCalled")
  private List<String> m_hapsCalled = new ArrayList<>();
  @Expose
  @SerializedName("dips")
  private List<String> m_dips = new ArrayList<>();
  @Expose
  @SerializedName("hapsAvailable")
  private List<String> m_hapsAvailable = new ArrayList<>();


  @Nonnull
  public List<String> getGenes() {
    return m_genes;
  }

  public void setGenes(List<String> genes) {
    m_genes.addAll(genes);
  }

  @Nonnull
  public List<String> getHapsCalled() {
    return m_hapsCalled;
  }

  public void setHapsCalled(List<String> haps) {
    m_hapsCalled.addAll(haps);
  }

  @Nonnull
  public List<String> getHapsAvailable() {
    return m_hapsAvailable;
  }

  public void setHapsAvailable(List<String> haps) {
    m_hapsAvailable.addAll(haps);
  }

  @Nonnull
  public List<String> getDips() {
    return m_dips;
  }

  public void setDips(List<String> dips) {
    m_dips.addAll(dips);
  }
}
