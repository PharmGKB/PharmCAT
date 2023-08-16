package org.pharmgkb.pharmcat.reporter.format.hl7;

public enum InteractionType {
  METABOLIZER("53040-2", "Genetic Variation's Effect on Drug Metabolism"),
  TRANSPORT("51961-1", "Genetic Variation's Effect on Drug Transport"),
  RISK("83009-1", "Genetic Variation's Effect on High-Risk"),
  EFFICACY("51961-1", "Genetic Variation's Effect on Drug Efficacy");

  private final String code;
  private final String header;

  InteractionType(String code, String header) {
    this.code = code;
    this.header = header;
  }

  public String getCode() {
    return this.code;
  }

  public String getHeader() {
    return this.header;
  }
}
