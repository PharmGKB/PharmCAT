
package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.List;
import javax.annotation.Nonnull;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;

public class Group implements Comparable<Group> {

  @Expose
  @SerializedName("id")
  private String id;
  @Expose
  @SerializedName("name")
  private String name;
  @Expose
  @SerializedName("annotations")
  private List<Annotation> annotations = new ArrayList<>();
  @Expose
  @SerializedName("genotypes")
  private List<String> genotypes = new ArrayList<>();
  @Expose
  @SerializedName("strength")
  private Strength strength;


  public String getId() {
    return id;
  }

  public void setId(String id) {
    this.id = id;
  }


  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }


  public List<Annotation> getAnnotations() {
    return annotations;
  }

  public void setAnnotations(List<Annotation> annotations) {
    this.annotations = annotations;
  }


  public List<String> getGenotypes() {
    return genotypes;
  }

  public void setGenotypes(List<String> genotypes) {
    this.genotypes = genotypes;
  }


  public Strength getStrength() {
    return strength;
  }

  public void setStrength(Strength strength) {
    this.strength = strength;
  }


  @Override
  public int compareTo(@Nonnull Group o) {

    if (id == null) {
      return -1;
    }
    else if (o.id == null) {
      return 1;
    }
    else {
      return id.compareTo(o.getId());
    }
  }
}
