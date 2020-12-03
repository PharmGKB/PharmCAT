
package org.pharmgkb.pharmcat.reporter.model;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;

/**
 * @deprecated
 */
public class RelatedGene {
  @Expose
  @SerializedName("id")
  private String id;
  @Expose
  @SerializedName("symbol")
  private String symbol;
  @Expose
  @SerializedName("name")
  private String name;


  public String getId() {
    return id;
  }

  public void setId(String id) {
    this.id = id;
  }

  public String getSymbol() {
    return symbol;
  }

  public void setSymbol(String symbol) {
    this.symbol = symbol;
  }

  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }
}
