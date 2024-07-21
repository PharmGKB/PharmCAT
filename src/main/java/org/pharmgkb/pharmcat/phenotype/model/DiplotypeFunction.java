package org.pharmgkb.pharmcat.phenotype.model;

import java.util.Map;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class DiplotypeFunction {
  @SerializedName("name")
  @Expose
  String name;
  @SerializedName("lookupKey")
  @Expose
  Map<String, Integer> lookupKey;
  @SerializedName("phenotype")
  @Expose
  String phenotype;
  @SerializedName("activityScore")
  @Expose
  String activityScore;


  public String getName() {
    return name;
  }

  public Map<String, Integer> getLookupKey() {
    return lookupKey;
  }

  public String getPhenotype() {
    return phenotype;
  }

  public String getActivityScore() {
    return activityScore;
  }
}
