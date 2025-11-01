package org.pharmgkb.pharmcat.reporter;

import java.util.List;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.SampleAllele;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;

import static org.junit.jupiter.api.Assertions.*;


/**
 * This is a JUnit test for {@link VariantReportFactory}.
 *
 * @author Mark Woon
 */
class VariantReportFactoryTest {


  @Test
  void renameSuballeles() throws Exception {

    Env env = new Env();

    DefinitionFile definitionFile = env.getDefinitionReader().getDefinitionFile("SLCO1B1");
    NamedAllele star45_001 = definitionFile.getNamedAllele("*45.001");
    assertNotNull(star45_001);
    VariantLocus vl = definitionFile.getVariantForPosition(star45_001.getCorePositions().first());

    VariantReportFactory vrFactory = new VariantReportFactory("SLCO1B1", "chr12", env);
    VariantReport vr = vrFactory.make(vl);
    assertTrue(vr.getAlleles().contains("*45"));
    assertFalse(vr.getAlleles().contains("*45.001"));


    SampleAllele sampleAllele = new SampleAllele("chr12", 21178615, "T", "C", true,
        true, null, List.of("T", "C"), "0|1", null, false);
    Variant variant = new Variant(vl, sampleAllele);

    vr = vrFactory.make(variant);
    assertTrue(vr.getAlleles().contains("*45"));
    assertFalse(vr.getAlleles().contains("*45.001"));
  }
}
