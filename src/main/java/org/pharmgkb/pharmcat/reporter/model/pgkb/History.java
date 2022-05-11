package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.Date;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * History event object for tracking history on pharmgkb objects
 *
 * @author Ryan Whaley
 */
public class History {

  @Expose
  @SerializedName("id")
  private long m_id;
  @Expose
  @SerializedName("date")
  private Date m_date;
  @Expose
  @SerializedName("type")
  private String m_type;
  @Expose
  @SerializedName("description")
  private String m_description;

  public long getId() {
    return m_id;
  }

  public Date getDate() {
    return m_date;
  }

  public String getType() {
    return m_type;
  }

  public String getDescription() {
    return m_description;
  }
}
