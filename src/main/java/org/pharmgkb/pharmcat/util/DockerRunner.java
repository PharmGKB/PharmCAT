package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.nio.file.attribute.PosixFilePermission;
import java.nio.file.attribute.PosixFilePermissions;
import java.util.Set;
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
   * Gets a temp mount directory for Docker to use.
   * If calling on a POSIX filesystem, this will make sure that the directory's permissions are set to 777 to  allow
   * Docker to work with it.
   */
  public static Path getTempMountDir() throws IOException {
    Path tmpDir = Files.createTempDirectory("pcat_docker-");
    tmpDir.toFile().deleteOnExit();
    if (FileSystems.getDefault().supportedFileAttributeViews().contains("posix")) {
      try {
        // make sure dir is writable
        Set<PosixFilePermission> perms = PosixFilePermissions.fromString("rwxrwxrwx");
        Files.setPosixFilePermissions(tmpDir, perms);
      } catch (Exception ex) {
        throw new IOException("Could not set permissions for " + tmpDir, ex);
      }
    }
    return tmpDir;
  }


  /**
   * Runs {@code bcftools norm} (multiallelic) via Docker.
   * Expects a container tagged "pcat" to be available.
   */
  public static void normalizeVcf(Path inFile, Path outFile) throws IOException {
    String dockerCmd = getDockerCmd(inFile.getParent());
    String toolCmd = "bcftools norm --no-version -m+ -c ws -Oz -f reference.fna.bgz -o data/" +
        PathUtils.getFilename(outFile) + " data/" + PathUtils.getFilename(inFile);

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
    // rename the copy back to the original file
    Files.move(bakFile, inFile);
    return bgzFile;
  }


  /**
   * Runs {@code pcat.prep_pharmcat_positions} via Docker.
   * This will index the pharmcat_positions.vcf.bgz and generate the uniallelic positions.
   * This expects that a container tagged "pcat" is available.
   */
  public static void prepPharmcatPositions(Path inFile) throws IOException {
    String dockerCmd = getDockerCmd(inFile.getParent());
    String localPath = "data/" + PathUtils.getFilename(inFile);

    runCmd(dockerCmd, "python -c \"import pcat; from pathlib import Path; " +
        "file = Path('" + localPath + "'); pcat.prep_pharmcat_positions(file, force_update=True, verbose=1)\"",
        inFile.getParent(), true);
  }


  private static void runCmd(String dockerCmd, String toolCmd, Path workDir) throws IOException {
    runCmd(dockerCmd, toolCmd, workDir, false);
  }

  private static void runCmd(String dockerCmd, String toolCmd, Path workDir, boolean verbose) throws IOException {
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
      if (!stdin.isEmpty()) {
        if (verbose) {
          System.out.println(stdin);
        } else {
          sf_logger.info(stdin);
        }
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
