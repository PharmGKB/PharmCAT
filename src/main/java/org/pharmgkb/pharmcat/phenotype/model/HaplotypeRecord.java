package org.pharmgkb.pharmcat.phenotype.model;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.jspecify.annotations.Nullable;

import static org.pharmgkb.pharmcat.reporter.TextConstants.isUnspecified;


public class HaplotypeRecord {
  @Expose
  @SerializedName("name")
  private String name;
  @Expose
  @SerializedName("activityValue")
  private @Nullable String activityValue;
  @Expose
  @SerializedName("functionValue")
  private @Nullable String functionValue;
  @Expose
  @SerializedName("lookupKey")
  private String lookupKey;


  public HaplotypeRecord(String name, @Nullable String activityValue, @Nullable String functionValue,
      @Nullable String lookupKey) {
    this.name = name;
    this.activityValue = activityValue;
    this.functionValue = functionValue;
    this.lookupKey = lookupKey;
  }

  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }

  public @Nullable String getActivityValue() {
    return activityValue;
  }

  // needed to support subsetting
  public void setActivityValue(String activityValue) {
    this.activityValue = activityValue;
    this.updateLookupKey();
  }

  public @Nullable String getFunctionValue() {
    return functionValue;
  }

  // needed to support subsetting
  public void setFunctionValue(String functionValue) {
    this.functionValue = functionValue;
    this.updateLookupKey();
  }

  public String getLookupKey() {
    return lookupKey;
  }

  public void setLookupKey(String lookupKey) {
    this.lookupKey = lookupKey;
  }

  private void updateLookupKey() {
    this.lookupKey = this.activityValue == null ? this.functionValue : this.activityValue;
  }

  public String toFormattedFunction() {
    if (isUnspecified(this.activityValue)) {
      return this.functionValue;
    } else {
      return String.format("Activity Value %s (%s)", this.activityValue, this.functionValue);
    }
  }

  public String toString() {
    return this.name;
  }
}
