package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import com.google.common.base.Splitter;
import org.apache.commons.lang3.ObjectUtils;


/**
 * Object to hold information about varaints from the tsv file.
 *
 * @author nate
 */
public class Variant implements Comparable<Variant> {
  private static final Splitter sf_semicolonSplitter = Splitter.on(";").trimResults().omitEmptyStrings();
  private String m_geneName;
  private String m_chromosome;
  private int m_position;
  private String m_genePosition;
  private ArrayList<String> HGVSg = new ArrayList<>();
  private ArrayList<String> m_cDna = new ArrayList<>();
  private ArrayList<String> m_proteinEffect = new ArrayList<>();
  private ArrayList<String> ALTs = new ArrayList<>();
  private String REF;
  private String m_rsid;


  public Variant(String chr, String geneName) {
    m_chromosome = chr;
    m_geneName = geneName;
  }


  public String getGeneName() {
    return m_geneName;
  }

  public String getChromosome() {
    return m_chromosome;
  }

  public int getPosition() {
    return m_position;
  }

  public void setPosition(String pos) {
    m_position = Integer.parseInt(pos);
  }


  public void setCDna(String cdna) {
    if (cdna.contains(";")) {
      System.out.println("??");
    }
    sf_semicolonSplitter.split(cdna).forEach(m_cDna::add);
  }

  public void addProteinEffect(String proteinEffect) {
    sf_semicolonSplitter.split(proteinEffect).forEach(m_proteinEffect::add);
  }

  public void addHGVSg(String hgvs) {
    sf_semicolonSplitter.split(hgvs).forEach(HGVSg::add);
  }

  public String getHGVSg() {
    return HGVSg.toString();
  }



  public String getGenePosition() {
    return m_genePosition;
  }

  public void setGenePosition(String genePosition) {
    m_genePosition = genePosition;
  }




  private int getStartPOS(String _HGVSg) {

    Pattern p = Pattern.compile("(\\d+)");
    Matcher m = p.matcher(_HGVSg);
    //System.out.println(_HGVSg);
    if (m.find()) {
      return Integer.parseInt(m.group(1));
    } else {
      return -1;
    }

  }

  public boolean setStartPOS() {
    if (HGVSg.isEmpty()) {
      return false;
    } else {
      m_position = getStartPOS(HGVSg.get(0));
      return true;
    }
  }

  public String getREF() {
    return REF;
  }

  public void setREF(String _REF) {
    REF = (_REF);
  }

  public String getRsid() {
    return m_rsid;
  }

  public void setRsid(String rsid) {
    m_rsid = rsid;
  }

  public void addALT(String _ALT) {
    ALTs.add(_ALT);
  }

  public ArrayList<String> getALTs() {
    return ALTs;
  }


  @Override
  public int compareTo(@Nonnull Variant o) {

    int rez = ObjectUtils.compare(m_chromosome, o.m_chromosome);
    if (rez != 0) {
      return rez;
    }
    rez = ObjectUtils.compare(m_position, o.m_position);
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_geneName, o.m_geneName);
  }
}
