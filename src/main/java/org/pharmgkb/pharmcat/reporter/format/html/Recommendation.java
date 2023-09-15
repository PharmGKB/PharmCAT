package org.pharmgkb.pharmcat.reporter.format.html;

import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import com.google.common.base.Preconditions;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.caller.DpydCaller;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.cpic.Publication;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationReport;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;


/**
 * This class collects all the {@link Report}s for a single drug.
 *
 * @author Mark Woon
 */
public class Recommendation implements Comparable<Recommendation> {
  private final Env m_env;
  private final String m_drug;
  private Report m_cpicReport;
  private Report m_dpwgReport;
  private boolean m_hasInferred;
  private boolean m_hasDpydInferred;
  private final Set<MessageAnnotation> m_messages = new LinkedHashSet<>();
  private final Set<MessageAnnotation> m_footnotes = new LinkedHashSet<>();
  private final List<Publication> m_citations = new ArrayList<>();


  public Recommendation(Env env, String drug) {
    m_env = env;
    m_drug = drug;
  }


  public void addReport(DataSource source, DrugReport report) {
    Preconditions.checkArgument(m_drug.equals(report.getName()));

    switch (source) {
      case CPIC -> {
        if (m_cpicReport != null) {
          throw new IllegalStateException("Multiple drug reports for " + report.getName() + " from " +
              source.getDisplayName());
        }
        m_cpicReport = new Report(source, report);
      }
      case DPWG -> {
        if (m_dpwgReport != null) {
          throw new IllegalStateException("Multiple drug reports for " + report.getName() + " from " +
              source.getDisplayName());
        }
        m_dpwgReport = new Report(source, report);
      }
    }

    boolean callMultiMatch = false;
    boolean scoreMultimatch = false;
    Multimap<String, Diplotype> mismatchDiplotypes = HashMultimap.create();
    for (GuidelineReport guideline : report.getGuidelines()) {
      for (GeneReport gr : guideline.getGeneReports()) {
        // pull in DPYD HapB3 warnings
        if (gr.getGene().equals(DpydCaller.GENE)) {
          for (MessageAnnotation ma : gr.getMessages()) {
            if (MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC.equals(ma.getName()) ||
                MessageHelper.MSG_DPYD_HAPB3_EXONIC_ONLY.equals(ma.getName())) {
              m_messages.add(ma);
            }
          }
        }
      }
      for (AnnotationReport annotation : guideline.getAnnotations()) {
        if (annotation.getGenotypes().size() > 1) {
          callMultiMatch = true;
        }
        for (Genotype genotype : annotation.getGenotypes()) {
          for (Diplotype diplotype : genotype.getDiplotypes()) {
            if (!scoreMultimatch && diplotype.getLookupKeys().size() > 1) {
              GenePhenotype gp = m_env.getPhenotype(diplotype.getGene(), source);
              scoreMultimatch = gp != null && gp.isMatchedByActivityScore();
            }
            if (diplotype.getOutsidePhenotypeMismatch() != null ||
                diplotype.getOutsideActivityScoreMismatch() != null) {
              mismatchDiplotypes.put(diplotype.getGene(), diplotype);
            }
          }
          if (scoreMultimatch || genotype.getDiplotypes().stream().anyMatch((d) -> {
            if (d.getLookupKeys().size() > 1) {
              GenePhenotype gp = m_env.getPhenotype(d.getGene(), source);
              return gp != null && gp.isMatchedByActivityScore();
            }
            return false;
          })) {
            scoreMultimatch = true;
          }
          if (genotype.isInferred()) {
            if (genotype.getDiplotypes().stream()
                .map(Diplotype::getGene)
                .anyMatch(g -> g.equals("DPYD"))) {
              m_hasDpydInferred = true;
            } else {
              m_hasInferred = true;
            }
            break;
          }
        }
      }
    }
    if (callMultiMatch) {
      m_messages.add(m_env.getMessage(MessageHelper.MSG_MULTI_CALL));
    }
    if (scoreMultimatch) {
      m_messages.add(m_env.getMessage(MessageHelper.MSG_MULTI_SCORE));
    }
    for (String gene : mismatchDiplotypes.keySet()) {
      StringBuilder builder = new StringBuilder()
          .append("Conflicting outside call data was provided for ")
          .append(gene)
          .append(".  PharmCAT will use provided ");
      if (mismatchDiplotypes.get(gene).stream().anyMatch(Diplotype::hasActivityScore)) {
        builder.append("activity score");
      } else {
        builder.append("phenotype");
      }
      builder.append(" to match recommendations.");
      m_messages.add(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "warn.mismatch.outsideCall",
          builder.toString()));
    }

    // messages
    report.getMessages().stream()
        .filter(MessageAnnotation.isMessage)
        .forEach(m_messages::add);
    // footnotes
    report.getMessages().stream()
        .filter(MessageAnnotation.isFootnote)
        .forEach(m_footnotes::add);
    // citations
    if (report.getCitations() != null && !report.getCitations().isEmpty()) {
      m_citations.addAll(report.getCitations());
    }
  }


  public String getDrug() {
    return m_drug;
  }


  public boolean isMatched() {
    return (m_cpicReport != null && m_cpicReport.isMatched()) ||
        (m_dpwgReport != null && m_dpwgReport.isMatched());
  }


  public List<Report> getReports() {
    List<Report> reports = new ArrayList<>();
    if (m_cpicReport != null) {
      reports.add(m_cpicReport);
    }
    if (m_dpwgReport != null) {
      reports.add(m_dpwgReport);
    }
    return reports;
  }

  public Report getCpicReport() {
    return m_cpicReport;
  }

  public Report getDpwgReport() {
    return m_dpwgReport;
  }

  public boolean isHasInferred() {
    return m_hasInferred;
  }

  public boolean isHasDpydInferred() {
    return m_hasDpydInferred;
  }

  public Set<MessageAnnotation> getMessages() {
    return m_messages;
  }

  public Set<MessageAnnotation> getFootnotes() {
    return m_footnotes;
  }

  public List<Publication> getCitations() {
    return m_citations;
  }


  @Override
  public int compareTo(Recommendation o) {
    return m_drug.compareTo(o.getDrug());
  }
}
