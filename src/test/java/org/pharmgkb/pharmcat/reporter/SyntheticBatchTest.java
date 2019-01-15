package org.pharmgkb.pharmcat.reporter;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.Map;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.PharmCAT;
import org.pharmgkb.pharmcat.VcfTestUtils;


/**
 * This test generates a collection of synthetic sample VCF files from our test resources and then runs the matcher and
 * reporter on those samples, writing the output to a desired directory.
 *
 * This is especially useful to see how different allele calls will be called and displayed in the final output report.
 *
 * @author Ryan Whaley
 */
public class SyntheticBatchTest {
  private static final Path sf_astrolabe = PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/test.astrolabe.tsv");
  private static final Map<String,String[]> sf_testVcfs = new LinkedHashMap<>();
  static {

    sf_testVcfs.put("cftr.reg_inc", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("cftr.ref_inc", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refG542X.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("cftr.inc_inc", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XG542X.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("cftr.ref_ref", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("example", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s1.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s1s1.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("cftr.reg_reg", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/F508delF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("cftr.ref_reg", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });



    sf_testVcfs.put("slco1b1.17.21", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s17s21.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("slco1b1.5.15", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("slco1b1.missing", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("slco1b1.multi", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/multi.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("slco1b1.no_match", new String[]{
        "SLCO1B1/multi.vcf"
    });

    sf_testVcfs.put("slco1b1.match_1", new String[]{
        "SLCO1B1/s5s15.vcf"
    });

    sf_testVcfs.put("cyp2c19.onlyRs12769205", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/rs12769205only.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("cyp2c19.refRs12769205", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/rs12769205ref.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("cyp2c19.s2s3", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/refref.vcf",
        "CYP2C19/s2s3.vcf",
        "CYP2C9/s1s1.vcf",
        "SLCO1B1/s1as1a.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("ugt1a1.phased.multi", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s60s80phased.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/F508delF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("ugt1a1.unphased.multi", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s60s80unphased.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/F508delF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("missing.genes", new String[]{
        "DPYD/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/F508delF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("cyp2c19.nonnormal", new String[]{
        "DPYD/s1s1.vcf",
        "TPMT/s1s1.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/F508delF508del.vcf",
        "CYP2C19/s1s35.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("tpmt.star1s", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1ss1ss3.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("tpmt.s1s1s", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1s1s.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("tpmt.hom1s_het3a", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/hom1s_het3a.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("tpmt.het3a", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/het3a.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("tpmt.s15offdata", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s15offdata.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("tpmt.unexpected.allele", new String[]{
        "TPMT/unexpected.allele.vcf"
    });

    sf_testVcfs.put("tpmt.unexpected.allele.with.s9", new String[]{
        "TPMT/unexpected.allele.with.s9.vcf"
    });

    sf_testVcfs.put("cyp2c19.rs28399504missing", new String[]{
        "DPYD/s1s1.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1ss1ss3.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s4bs17rs28399504missing.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });


    sf_testVcfs.put("dpyd.stars12b", new String[]{
        "DPYD/s1s2b.vcf",
        "UGT1A1/s1s1.vcf",
        "TPMT/s1ss1ss3.vcf",
        "CYP3A5/s1s7.vcf",
        "CFTR/G542XF508del.vcf",
        "CYP2C19/s2s2.vcf",
        "CYP2C9/s2s3.vcf",
        "SLCO1B1/s5s15.vcf",
        "VKORC1/-1639A-1639A.vcf",
        "cyp4f2/s1s1.vcf",
        "IFNL3/rs12979860CC.vcf"
    });

    sf_testVcfs.put("cyp2c19.rxPossible", new String[]{
        "CYP2C19/s1s17.vcf"
    });


    sf_testVcfs.put("ugt1a1.s1s1", new String[]{
      "UGT1A1/s1s1.vcf"
    });

    sf_testVcfs.put("ugt1a1.s1s28s60s80unphased", new String[]{
      "UGT1A1/s1s28s60s80unphased.vcf"
    });

    sf_testVcfs.put("ugt1a1.s1s60s80unphased", new String[]{
      "UGT1A1/s1s60s80unphased.vcf"
    });

    sf_testVcfs.put("ugt1a1.s1s60s80phased", new String[]{
      "UGT1A1/s1s60s80phased.vcf"
    });

    sf_testVcfs.put("ugt1a1.s6s6", new String[]{
      "UGT1A1/s6s6.vcf"
    });

    sf_testVcfs.put("ugt1a1.s28s37", new String[]{
      "UGT1A1/s28s37.vcf"
    });

    sf_testVcfs.put("ugt1a1.s28s80unphased", new String[]{
      "UGT1A1/s28s80unphased.vcf"
    });

    sf_testVcfs.put("ugt1a1.s28s80phased", new String[]{
      "UGT1A1/s28s80phased.vcf"
    });

    sf_testVcfs.put("ugt1a1.phased.s28s80s6s60", new String[]{
        "UGT1A1/s28s80s6s60phased.vcf"
    });

    sf_testVcfs.put("ugt1a1.unphased.s28s80s6s60", new String[]{
        "UGT1A1/s28s80s6s60unphased.vcf"
    });

    sf_testVcfs.put("ugt1a1.s6s60s80s28missingphased", new String[]{
        "UGT1A1/s6s60s80s28missingphased.vcf"
    });

    sf_testVcfs.put("ugt1a1.s6s60s80s28missingunphased", new String[]{
        "UGT1A1/s6s60s80s28missingunphased.vcf"
    });

    sf_testVcfs.put("ugt1a1.s80s28missing", new String[]{
        "UGT1A1/s80s28missing.vcf"
    });

    sf_testVcfs.put("ugt1a1.s28het", new String[]{
        "UGT1A1/s28het.vcf"
    });

    sf_testVcfs.put("ugt1a1.na12717", new String[]{
        "UGT1A1/NA12717_UGT1A1.vcf"
    });

    sf_testVcfs.put("ugt1a1.na18868", new String[]{
        "UGT1A1/NA18868_UGT1A1.vcf"
    });

    sf_testVcfs.put("ugt1a1.na19785", new String[]{
        "UGT1A1/NA19785_UGT1A1.vcf"
    });

    sf_testVcfs.put("ugt1a1.s28s28unphaseds60s80miss", new String[]{
        "UGT1A1/s28s28unphaseds60s80miss.vcf"
    });

    sf_testVcfs.put("ugt1a1.s27s28unphaseds80s60missing", new String[]{
        "UGT1A1/s27s28unphaseds80s60missing.vcf"
    });

    sf_testVcfs.put("ugt1a1.s28hets60homounphaseds80missing", new String[]{
        "UGT1A1/s28hets60homounphaseds80missing.vcf"
    });

    sf_testVcfs.put("ugt1a1.HG00436", new String[]{
        "UGT1A1/HG00436.vcf"
    });

    sf_testVcfs.put("ugt1a1.s1s80s27s60s28missingphased", new String[]{
        "UGT1A1/s1s80s27s60s28missingphased.vcf"
    });

    sf_testVcfs.put("ugt1a1.s1s60s80s6phased", new String[]{
        "UGT1A1/s1s60s80s6phased.vcf"
    });

    sf_testVcfs.put("ugt1a1.s1s60s80s28s6phased", new String[]{
        "UGT1A1/s1s60s80s28s6phased.vcf"
    });

    sf_testVcfs.put("ugt1a1.s1s37s80s60phased", new String[]{
        "UGT1A1/s1s37s80s60phased.vcf"
    });


    sf_testVcfs.put("cyp2c19.rs12248560missing", new String[]{
        "cyp2c19/s1s1rs12248560missing.vcf"
    });

    sf_testVcfs.put("cyp3a5.s1s1rs776746missing", new String[]{
        "cyp3a5/s1s1rs776746missing.vcf"
    });

    sf_testVcfs.put("cftr.G1244Eref", new String[]{
        "cftr/G1244Eref.vcf"
    });

    sf_testVcfs.put("cftr.G1244EF508del", new String[]{
        "cftr/G1244EF508del.vcf"
    });

    sf_testVcfs.put("cftr.G551DG542X", new String[]{
        "cftr/G551DG542X.vcf"
    });

    sf_testVcfs.put("cftr.inc_vwithrecom", new String[]{
        "cftr/inc_vwithrecom.vcf"
    });

    sf_testVcfs.put("cftr.ref_vwithrecom", new String[]{
        "cftr/ref_vwithrecom.vcf"
    });

    sf_testVcfs.put("cftr.v_vwithrecom", new String[]{
        "cftr/v_vwithrecom.vcf"
    });

    sf_testVcfs.put("cftr.v_vwithrecom1", new String[]{
        "cftr/v_vwithrecom1.vcf"
    });

    sf_testVcfs.put("cyp2c9.s2s24", new String[]{
        "cyp2c9/s2s24.vcf"
    });

    sf_testVcfs.put("cftr.refI507missing", new String[]{
        "cftr/refI507missing.vcf"
    });

    sf_testVcfs.put("cftr.F508delfirstsecondmis", new String[]{
        "cftr/F508delfirstsecondmis.vcf"
    });

    sf_testVcfs.put("cftr.F508delfirstsecondref", new String[]{
        "cftr/F508delfirstsecondref.vcf"
    });

    sf_testVcfs.put("cftr.F508delsecondreffirst", new String[]{
        "cftr/F508delsecondreffirst.vcf"
    });

    sf_testVcfs.put("cftr.F508delsecondfirstmissing", new String[]{
        "cftr/F508delsecondfirstmissing.vcf"
    });

    sf_testVcfs.put("cyp2c9.s2s24only", new String[]{
        "cyp2c9/s2s24only.vcf"
    });

    sf_testVcfs.put("cyp2c19.rs12769205call", new String[]{
        "cyp2c19/rs12769205call.vcf"
    });

    sf_testVcfs.put("cyp2c19.s4s17het", new String[]{
        "cyp2c19/s4s17het.vcf"
    });

    sf_testVcfs.put("dpyd.c557c703", new String[]{
        "DPYD/c557c703.vcf"
    });

    sf_testVcfs.put("dpyd.c1129-5923c2846", new String[]{
        "DPYD/c1129-5923c2846.vcf"
    });

    sf_testVcfs.put("dpyd.c1156het", new String[]{
        "DPYD/c1156het.vcf"
    });

    sf_testVcfs.put("dpyd.c1679c1156", new String[]{
        "DPYD/c1679c1156.vcf"
    });

    sf_testVcfs.put("dpyd.c2846het", new String[]{
        "DPYD/c2846het.vcf"
    });

    sf_testVcfs.put("dpyd.novariant", new String[]{
        "DPYD/novariant.vcf"
    });
  }

  private PharmCAT m_pharmcat;
  private Path m_outputDir;

  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("o", "output-dir", "directory to output to", true, "o")
        .addOption("g", "guideline-dir", "directory of guideline annotations (JSON files)", false, "n");

    try {
      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path outputDir = cliHelper.getValidDirectory("o", true);
      Path guidelineDir = null;
      if (cliHelper.hasOption("g")) {
        guidelineDir = cliHelper.getValidDirectory("g", false);
      }

      SyntheticBatchTest piplelineTest = new SyntheticBatchTest(outputDir, guidelineDir);
      piplelineTest.execute();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }


  private SyntheticBatchTest(Path outputDir, Path guidelineDir) throws IOException {
    m_outputDir = outputDir;
    m_pharmcat = new PharmCAT(outputDir, null, guidelineDir).keepMatcherOutput();
    m_pharmcat.writeJson(true);
  }

  private void execute() throws Exception {
    System.out.println("Run time: " + new Date());

    for (String key : sf_testVcfs.keySet()) {
      Path sampleDir = m_outputDir.resolve(key);
      if (!sampleDir.toFile().exists()) {
        if (!sampleDir.toFile().mkdirs()) {
          throw new RuntimeException("Output directory could not be created " + sampleDir.toAbsolutePath());
        }
      }
      
      Path sampleVcf = writeVcf(sampleDir.resolve(key + ".vcf"), sf_testVcfs.get(key));
      m_pharmcat.setOutputDir(sampleDir);
      m_pharmcat.execute(sampleVcf, sf_astrolabe, null);
      System.out.println("Generated " + key);
    }
  }

  private Path writeVcf(Path outputVcf, String[] filesToInclude) {
    try (FileWriter writer = new FileWriter(outputVcf.toFile())) {
      writer.write(VcfTestUtils.writeVcf(filesToInclude));
    } catch (IOException e) {
      e.printStackTrace();
    }

    return outputVcf;
  }
}
