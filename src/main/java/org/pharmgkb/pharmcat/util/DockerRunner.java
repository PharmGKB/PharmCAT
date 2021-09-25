package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.nio.file.Path;
import org.pharmgkb.common.util.PathUtils;


/**
 * Helper class to run docker.
 *
 * @author Mark Woon
 */
public class DockerRunner {
  private static final boolean sf_isWindows = System.getProperty("os.name").toLowerCase().startsWith("windows");


  /**
   * Runs {@code bcftools norm} via Docker.  Expects a container tagged "pcat" to be available.
   */
  public static void normalizeVcf(Path inFile, Path outFile) throws IOException {
    String dockerCmd = "docker run --rm -v " + inFile.getParent() + ":/pharmcat/data pcat ";
    String bcfToolsCmd = "bcftools norm -m+ -c ws -Oz -f reference.fasta -o data/" + PathUtils.getFilename(outFile) +
        " data/" + PathUtils.getFilename(inFile);

    ProcessBuilder builder = new ProcessBuilder();
    if (sf_isWindows) {
      builder.command("cmd.exe", "/c", dockerCmd + bcfToolsCmd);
    } else {
      builder.command("sh", "-c", dockerCmd + bcfToolsCmd);
    }
    builder.directory(inFile.getParent().toFile());
    Process process = builder.start();
    try {
      int exitCode = process.waitFor();
      String stdin = new String(process.getInputStream().readAllBytes());
      if (stdin.length() > 0) {
        System.out.println("OUTPUT:");
        System.out.println(stdin);
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
