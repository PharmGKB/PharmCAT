package org.pharmgkb.pharmcat;

import org.checkerframework.checker.nullness.qual.Nullable;


/**
 * This is an enumeration of possible {@link Pipeline} results.
 *
 * @author Mark Woon
 */
public class PipelineResult {
  public enum Status {
    NOOP,
    SUCCESS,
    FAILURE
  }
  private final Status m_status;
  private final String m_sampleId;
  private final String m_basename;


  public PipelineResult(Status status, String basename, @Nullable String sampleId) {
    m_status = status;
    m_basename = basename;
    m_sampleId = sampleId;
  }


  public Status getStatus() {
    return m_status;
  }


  public String getBasename() {
    return m_basename;
  }


  public @Nullable String getSampleId() {
    return m_sampleId;
  }
}
