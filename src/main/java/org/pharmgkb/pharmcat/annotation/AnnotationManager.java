package org.pharmgkb.pharmcat.annotation;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Properties;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.SheetsHelper;


/**
 * This class manages annotation files.
 *
 * @author Mark Woon
 */
public class AnnotationManager {
  private String m_googleUser;
  private String m_googleKey;


  private AnnotationManager(Path propertyFile) throws IOException {

    Properties properties = new Properties();
    try (BufferedReader reader = Files.newBufferedReader(propertyFile)) {
      properties.load(reader);
    }
    m_googleUser = StringUtils.stripToNull((String)properties.get("google.user"));
    Preconditions.checkState(m_googleUser != null, "Missing property: 'google.user");
    m_googleKey = StringUtils.stripToNull((String)properties.get("google.key"));
    Preconditions.checkState(m_googleKey != null, "Missing property: 'google.key");
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("p", "properties-file", "PharmCAT properties file", false, "p")
          .addOption("d", "download", "download annotation files", true, "d");

      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path propsFile = CliUtils.getPropsFile(cliHelper, "p");
      Path downloadDir = cliHelper.getValidDirectory("d", true);

      AnnotationManager manager = new AnnotationManager(propsFile);
      manager.download(downloadDir);

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  private void download(Path downloadDir) throws Exception {

    SheetsHelper sh = new SheetsHelper(m_googleUser, m_googleKey);
    sh.downloadRsidAnnotations(downloadDir.resolve(AnnotationReader.RSID_ANNOTATIONS_FILENAME));
    sh.downloadNamedAlleleAnnotations(downloadDir.resolve(AnnotationReader.NAMED_ALLELE_ANNOTATIONS_FILENAME));
  }
}
