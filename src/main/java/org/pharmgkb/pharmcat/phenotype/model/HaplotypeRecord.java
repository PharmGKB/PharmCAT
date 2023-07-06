package org.pharmgkb.pharmcat.phenotype.model;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;

import static org.pharmgkb.pharmcat.reporter.TextConstants.isUnspecified;


public class HaplotypeRecord {
  @Expose
  @SerializedName("name")
  private String name;
  @Expose
  @SerializedName("activityValue")
  private String activityValue;
  @Expose
  @SerializedName("functionValue")
  private String functionValue;
  @Expose
  @SerializedName("lookupKey")
  private String lookupKey;

  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }

  public String getActivityValue() {
    return activityValue;
  }

  public void setActivityValue(String activityValue) {
    this.activityValue = activityValue;
  }

  public String getFunctionValue() {
    return functionValue;
  }

  public void setFunctionValue(String functionValue) {
    this.functionValue = functionValue;
  }

  public String getLookupKey() {
    return lookupKey;
  }

  public void setLookupKey(String lookupKey) {
    this.lookupKey = lookupKey;
  }

  public String toFormattedFunction() {
    if (isUnspecified(this.activityValue)) {
      return this.functionValue;
    } else {
      return String.format("%s (%s)", this.activityValue, this.functionValue);
    }
  }

  public String toString() {
    return this.name;
  }
}
