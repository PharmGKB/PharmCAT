package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Set;
import java.util.jar.JarFile;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class DistributionArchiveTest {

  @Test
  void thinArchiveKeepsExistingContents() throws IOException {
    Set<String> entries = readEntries("pharmcat.thinJar");

    assertTrue(entries.contains("org/pharmgkb/pharmcat/PharmCAT.class"));
    assertFalse(entries.stream().anyMatch(name -> name.startsWith("ch/qos/logback/")));
    assertFalse(entries.stream().anyMatch(name -> name.startsWith("org/slf4j/")));
    assertTrue(entries.contains("logback.xml"));
  }

  @Test
  void executableArchiveKeepsBundledLoggingBehavior() throws IOException {
    Set<String> entries = readEntries("pharmcat.shadowJar");

    assertTrue(entries.contains("org/pharmgkb/pharmcat/PharmCAT.class"));
    assertTrue(entries.contains("ch/qos/logback/classic/Logger.class"));
    assertTrue(entries.contains("logback.xml"));
  }

  @Test
  void noLogbackArchiveBundlesDependenciesWithoutLogback() throws IOException {
    String property = "pharmcat.noLogbackJar";
    Set<String> entries = readEntries(property);

    assertTrue(archivePath(property).getFileName().toString().endsWith("-no-logback.jar"));
    assertTrue(entries.contains("org/pharmgkb/pharmcat/PharmCAT.class"));
    assertTrue(entries.contains("com/google/gson/Gson.class"));
    assertTrue(entries.contains("org/slf4j/Logger.class"));
    assertFalse(entries.stream().anyMatch(name -> name.startsWith("ch/qos/logback/")));
    assertFalse(entries.contains("logback.xml"));

    try (JarFile jarFile = new JarFile(archivePath(property).toFile())) {
      assertEquals("org.pharmgkb.pharmcat.PharmCAT",
          jarFile.getManifest().getMainAttributes().getValue("Main-Class"));
    }
  }

  private static Set<String> readEntries(String property) throws IOException {
    Path archive = archivePath(property);
    try (JarFile jarFile = new JarFile(archive.toFile())) {
      return jarFile.stream().map(entry -> entry.getName()).collect(Collectors.toSet());
    }
  }

  private static Path archivePath(String property) {
    return Path.of(System.getProperty(property));
  }
}
