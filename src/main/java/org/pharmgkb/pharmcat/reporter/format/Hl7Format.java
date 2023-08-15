package org.pharmgkb.pharmcat.reporter.format;

import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.StringJoiner;
import java.lang.StringBuilder;
import java.util.TreeSet;
import ca.uhn.hl7v2.DefaultHapiContext;
import ca.uhn.hl7v2.HL7Exception;
import ca.uhn.hl7v2.HapiContext;
import ca.uhn.hl7v2.model.Message;
import ca.uhn.hl7v2.parser.Parser;
import ca.uhn.hl7v2.util.Terser;

import com.google.common.collect.ImmutableMap;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationReport;
import org.pharmgkb.pharmcat.reporter.model.result.DrugLink;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * An HL7/FHIR formatted version of {@link ReportContext} data.
 */
public class Hl7Format extends AbstractFormat {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final ImmutableMap<String, String> sf_geneToInteractionType =
      new ImmutableMap.Builder<String, String>()
          .put("CYP2B6", "Metabolizer")
          .put("CYP2C9", "Metabolizer")
          .put("CYP2C19", "Metabolizer")
          .put("CYP2D6", "Metabolizer")
          .put("CYP3A5", "Metabolizer")
          .put("DPYD", "Metabolizer")
          .put("NUDT15", "Metabolizer")
          .put("TPMT", "Metabolizer")
          .put("UGT1A1", "Metabolizer")
          .put("ABCG2", "Transport")
          .put("SLCO1B1", "Transport")
          .put("CACNA1S", "Risk")
          .put("G6PD", "Risk")
          .put("HLA-A", "Risk")
          .put("HLA-B", "Risk")
          .put("MT-RNR1", "Risk")
          .put("RYR1", "Risk")
          .put("CFTR", "Efficacy")
          .build();
  private static final ImmutableMap<String, String> sf_interactionTypeToCPICCode =
      new ImmutableMap.Builder<String, String>()
          .put("Metabolizer", "53040-2")
          .put("Transport", "51961-1")
          .put("Risk", "83009-1")
          .put("Efficacy", "51961-1")
          .build();
  private static final ImmutableMap<String, String> sf_interactionTypeToCPICHeader =
      new ImmutableMap.Builder<String, String>()
          .put("Metabolizer", "Genetic Variation's Effect on Drug Metabolism")
          .put("Transport", "Genetic Variation's Effect on Drug Transport")
          .put("Risk", "Genetic Variation's Effect on High-Risk")
          .put("Efficacy", "Genetic Variation's Effect on Drug Efficacy")
          .build();
  private static final ImmutableMap<String, String> sf_susceptibilityToRiskLevel =
      new ImmutableMap.Builder<String, String>()
          .put("Malignant Hyperthermia Susceptibility", "High Risk")
          .put("Uncertain Susceptibility", "Normal Risk")
          .build();
  private static final ImmutableMap<String, String> sf_cftrGeneToEfficacy =  new ImmutableMap.Builder<String, String>()
      .put("ivacaftor responsive in CF patients", "Responsive")
      .put("ivacaftor non-responsive in CF patients", "Presumed non-responsive")
      .build();

  private final SimpleDateFormat m_dateFormat = new SimpleDateFormat("yyyyMMddHHmmss");
  private final String f_order;


  public Hl7Format(Path outputPath, Env env, Path orderPath) throws IOException {
    super(outputPath, env);
    f_order = Files.readString(orderPath);
  }


  public void write(ReportContext reportContext) throws IOException {

    try (BufferedWriter writer = Files.newBufferedWriter(getOutputPath(), StandardCharsets.UTF_8)) {
      writer.write(generateHL7(reportContext));
    }
  }

  private String generateHL7(ReportContext reportContext) {
    try (HapiContext context = new DefaultHapiContext()) {
      Parser parser = context.getGenericParser();
      Message message = parser.parse(f_order);
      Terser terser = new Terser(message);

      String sep = terser.get("MSH-1");
      String subSep = terser.get("MSH-2").substring(0, 1);

      String msh = generateMSHSegment(terser, sep, subSep);
      String pid = generatePIDSegment(terser, sep, subSep);
      String obr = generateOBRSegment(terser, sep, subSep);
      String obx = generateOBXSegments(reportContext, sep, subSep);

      return msh + pid + obr + obx;
    } catch (HL7Exception | IOException ex) {
      sf_logger.error("Error generating HL-7 report", ex);
    }
    return "";
  }

  private String getCPICRiskMessage(String gene, String pheno) {

    if (gene.equals("CACNA1S") || gene.equals("RYR1")) {
      return sf_susceptibilityToRiskLevel.get(pheno);
    }
    else if (gene.equals("MT-RNR1")) {
      if (pheno.contains("increased risk")) {
        return "High Risk";
      }
      else {
        return "Normal Risk";
      }
    }
    else if (gene.equals("G6PD")) {
      if (pheno.contains("Deficient") || pheno.equals("Variable")) {
        return "High Risk";
      }
      else {
        return "Normal Risk";
      }
    }
    else if (gene.equals("HLA-A") || gene.equals("HLA-B")) {
      // mapping provided for these terms was less concrete and based on whether the term contains the words positive
      // or negative
      if (pheno.contains("positive")) {
        return "High Risk";
      }
      else {
        return "Normal Risk";
      }
    }
    return "";
  }

  private String generateMSHSegment(Terser t, String sep, String subsep) throws HL7Exception {
    StringJoiner msh_joiner = new StringJoiner(sep);

    msh_joiner.add("MSH");
    msh_joiner.add(t.get("MSH-2"));
    msh_joiner.add(t.get("MSH-6")); //MSH-3: sending application
    msh_joiner.add(t.get("MSH-4"));
    msh_joiner.add("").add("");
    msh_joiner.add(m_dateFormat.format(new Date())); //MSH-7: date
    msh_joiner.add("");
    msh_joiner.add("ORU" + subsep + "R01"); //MSH-9-1: msg type, and MSH-9-2: trigger event
    msh_joiner.add(Long.toString(System.currentTimeMillis())); //MSH-10: control ID (server time in millis)
    msh_joiner.add("P"); //MSH-11: processing ID (P = production)
    msh_joiner.add("2.3").add("").add("\n"); //MSH-12: HL7 version, MSH-13 and MSH-14 left empty
    return msh_joiner.toString();
  }

  private String generatePIDSegment(Terser t, String sep, String subsep) throws HL7Exception {
    StringJoiner pid_joiner = new StringJoiner(sep);
    pid_joiner.add("PID");
    pid_joiner.add(t.get("PID-1"));
    pid_joiner.add("");
    pid_joiner.add(t.get("PID-2-1") + subsep.repeat(4) + t.get("PID-2-4")); //PID-3-1 and PID-3-5
    pid_joiner.add("");
    pid_joiner.add(t.get("PID-5-1") + subsep + t.get("PID-5-2"));
    pid_joiner.add("");
    pid_joiner.add(t.get("PID-7")); //PID-7: DOB
    pid_joiner.add(t.get("PID-8") + "\n"); //PID-8: gender
    return pid_joiner.toString();
  }

  private String generateOBRSegment(Terser t, String sep, String subsep) throws HL7Exception {
    StringJoiner obr_joiner = new StringJoiner(sep);
    obr_joiner.add("OBR");
    obr_joiner.add(t.get("OBR-1"));
    obr_joiner.add(t.get("OBR-2"));
    obr_joiner.add(Long.toString(System.currentTimeMillis()));
    obr_joiner.add(t.get("OBR-4-1") + subsep + t.get("OBR-4-2"));
    obr_joiner.add("").add("");
    obr_joiner.add(m_dateFormat.format(new Date())); //Observation date & time
    obr_joiner.add("").add("").add("").add("").add("").add("").add("");
    obr_joiner.add(t.get("OBR-16-1") + subsep + t.get("OBR-16-2") + subsep + t.get("OBR-16-3") +
        subsep + t.get("OBR-16-4") + subsep.repeat(9) + t.get("OBR-16-9"));
    obr_joiner.add("").add("").add("").add("").add("");
    obr_joiner.add(m_dateFormat.format(new Date()));
    obr_joiner.add("").add("");
    obr_joiner.add("F"); //status to Final
    obr_joiner.add("").add("");
    obr_joiner.add(t.get("OBR-28-1") + subsep + t.get("OBR-28-2") + subsep + t.get("OBR-28-3") + subsep +
        subsep.repeat(6) + t.get("OBR-28-9") + subsep.repeat(4) + t.get("OBR-28-13") + "\n");
    return obr_joiner.toString();
  }

  private String generateOBXSegments(ReportContext reportContext, String sep, String subsep) {

    Map<String, String> drugToRecommendation = new HashMap<>();
    // only consider CPIC results
    for (DrugReport drugReport : new TreeSet<>(reportContext.getDrugReports().get(DataSource.CPIC).values())) {
      if (drugReport.getGuidelines().size() == 0) {
        continue;
      }
      for (GuidelineReport guidelineReport : drugReport.getGuidelines()) {
        for (AnnotationReport report : guidelineReport.getAnnotations()) {
          drugToRecommendation.put(drugReport.getName(), report.getDrugRecommendation());
        }
      }
    }

    StringBuilder obx = new StringBuilder();
    String postfix = (new StringBuilder()).append(sep.repeat(5)).append("F\r\n").toString();
    String subId = "";

    int currentOBXSegment = 1;

    // only considering CPIC results
    for (GeneReport report : reportContext.getGeneReports().get(DataSource.CPIC).values()) {
      String gene = report.getGene();

      // TODO: you are assuming that there is only 1 diplotype, which is not correct
      // TODO: you will also need to distinguish between source and recommendation diplotypes
      String diplo = report.getRecommendationDiplotypes().first().toString();
      String pheno = String.join(", ", report.getRecommendationDiplotypes().first().getPhenotypes());

      String interactionType = sf_geneToInteractionType.get(gene);
      String code = sf_interactionTypeToCPICCode.get(interactionType);
      String header = sf_interactionTypeToCPICHeader.get(interactionType);

      String interactionMessage;
      if ("Risk".equals(interactionType)) {
        interactionMessage = getCPICRiskMessage(gene, pheno);
      }
      else if ("Efficacy".equals(interactionType)){
        interactionMessage = sf_cftrGeneToEfficacy.get(pheno);
      }
      else {
        interactionMessage = pheno;
      }
      obx.append("OBX").append(sep).append(currentOBXSegment++).append(sep).append("ST").append("47998-0^Variant Display Name^LN").append(sep).append(subId).append(sep).append(gene).append(postfix);
      obx.append("OBX").append(sep).append(currentOBXSegment++).append(sep).append("CWE").append("48018-6^Gene studied^LN").append(sep).append(subId).append(sep).append(gene).append(postfix);
      obx.append("OBX").append(sep).append(currentOBXSegment++).append(sep).append("ST").append("84413-4^Genotype display name^LN").append(sep).append(subId).append(sep).append(diplo).append(postfix);
      obx.append("OBX").append(sep).append(currentOBXSegment++).append(sep).append("ST").append(code).append(subsep).append(header).append(subsep).append("LN").append(sep).append(subId).append(sep).append(interactionMessage).append(postfix);

      for (DrugLink drugLink : report.getRelatedDrugs()) {
        obx.append("OBX").append(sep).append(currentOBXSegment++).append(sep).append("CWE").append("51963-7^Medication Assessed^LN").append(sep).append(subId).append(sep).append(drugLink.getName()).append(postfix);

        //How to populate medication usage type remains an open question
        obx.append("OBX").append(sep).append(currentOBXSegment++).append(sep).append("ST").append("82116-5^Medication Usage Suggestion [Type]^LN").append(sep).append(subId ).append(sep).append("type").append(postfix);
        obx.append("OBX").append(sep).append(currentOBXSegment++).append(sep).append("TX").append("83010-9^Medication Usage Suggestion [Narrative]^LN").append(sep).append(subId).append(sep).append(drugToRecommendation.get(drugLink.getName())).append(postfix);
      }
    }
    return obx.toString();
  }
}

