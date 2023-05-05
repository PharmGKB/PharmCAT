package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import org.pharmgkb.common.util.PathUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Helper class to run docker.
 * <p>
 * <b>Expects a container tagged "pcat" to be available.</b>
 *
 * @author Mark Woon
 */
public class DockerRunner {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final boolean sf_isWindows = System.getProperty("os.name").toLowerCase().startsWith("windows");


  private static String getDockerCmd(Path dataDir) {
    return "docker run --rm -v " + dataDir + ":/pharmcat/data pcat ";
  }


  /**
   * Runs {@code bcftools norm} via Docker.
   * Expects a container tagged "pcat" to be available.
   */
  public static void normalizeVcf(Path inFile, Path outFile) throws IOException {
    String dockerCmd = getDockerCmd(inFile.getParent());
    String toolCmd = "bcftools norm -m+ -c ws -Oz -f reference.fna.bgz -o data/" + PathUtils.getFilename(outFile) +
        " data/" + PathUtils.getFilename(inFile);

    runCmd(dockerCmd, toolCmd, inFile.getParent());
  }


  /**
   * Runs {@code bgzip} via Docker.
   * Expects a container tagged "pcat" to be available.
   */
  public static Path bgzip(Path inFile) throws IOException {
    String dockerCmd = getDockerCmd(inFile.getParent());
    // cannot get file redirect to work, so cannot use -c flag
    String toolCmd = "bgzip -f data/" + PathUtils.getFilename(inFile);

    // make a copy
    Path bakFile = inFile.getParent().resolve(PathUtils.getFilename(inFile) + ".bak");
    Files.copy(inFile, bakFile, StandardCopyOption.REPLACE_EXISTING);

    runCmd(dockerCmd, toolCmd, inFile.getParent());

    Path gzFile = inFile.getParent().resolve(PathUtils.getFilename(inFile) + ".gz");
    Path bgzFile = inFile.getParent().resolve(PathUtils.getFilename(inFile) + ".bgz");
    // rename .gz to .bgz
    Files.move(gzFile, bgzFile, StandardCopyOption.REPLACE_EXISTING);
    // rename copy to original file
    Files.move(bakFile, inFile);
    return bgzFile;
  }


  /**
   * Runs {@code bcftools index} via Docker.
   * Expects a container tagged "pcat" to be available.
   */
  public static void indexVcf(Path inFile) throws IOException {
    String dockerCmd = getDockerCmd(inFile.getParent());
    String localPath = "data/" + PathUtils.getFilename(inFile);

    runCmd(dockerCmd, "python -c \"import preprocessor; from pathlib import Path; " +
        "file = Path('" + localPath + "'); preprocessor.prep_pharmcat_positions(file)\"", inFile.getParent());
  }


  private static void runCmd(String dockerCmd, String toolCmd, Path workDir) throws IOException {

    ProcessBuilder builder = new ProcessBuilder();
    if (sf_isWindows) {
      builder.command("cmd.exe", "/c", dockerCmd + toolCmd);
    } else {
      builder.command("sh", "-c", dockerCmd + toolCmd);
    }
    builder.directory(workDir.toFile());
    Process process = builder.start();
    try {
      int exitCode = process.waitFor();
      String stdin = new String(process.getInputStream().readAllBytes());
      if (stdin.length() > 0) {
        sf_logger.info(stdin);
      }
      if (exitCode != 0) {
        throw new IOException("Docker run failed (" + exitCode + ")\n" +
            new String(process.getErrorStream().readAllBytes()));
      }
    } catch (InterruptedException ex) {
      Thread.currentThread().interrupt();
      throw new IOException("Docker run interrupted", ex);
    }
  }
}
