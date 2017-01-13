
package org.pharmgkb.pharmcat.reporter.model;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.model.result.Markdown;


public class Annotation {

  @SerializedName("id")
  @Expose
  private Integer id;
  @SerializedName("type")
  @Expose
  private Type type;
  @SerializedName("markdown")
  @Expose
  private Markdown m_markdown;

  /**
   *
   * @return
   *     The id
   */
  public Integer getId() {
    return id;
  }

  /**
   *
   * @param id
   *     The id
   */
  public void setId(Integer id) {
    this.id = id;
  }

  /**
   *
   * @return
   *     The type
   */
  public Type getType() {
    return type;
  }

  /**
   *
   * @param type
   *     The type
   */
  public void setType(Type type) {
    this.type = type;
  }

  public Markdown getMarkdown() {
    return m_markdown;
  }

  public void setMarkdown(Markdown markdown) {
    m_markdown = markdown;
  }
}
