package org.pharmgkb.pharmcat.util;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;


/**
 * Helper methods for working with zip files.
 *
 * @author Mark Woon
 */
public class ZipUtils {

  public static void unzip(Path zipFile, Path destDir) throws IOException {
    try (ZipInputStream zis = new ZipInputStream(new FileInputStream(zipFile.toFile()))) {
      // list files in zip
      ZipEntry zipEntry = zis.getNextEntry();
      while (zipEntry != null) {
        Path newFile = newFile(destDir, zipEntry);

        if (zipEntry.isDirectory()) {
          if (!Files.exists(newFile)) {
            Files.createDirectories(newFile);
          }
        } else {
          // some zips store file path only, need create parent directories
          Path parent = newFile.getParent();
          if (!Files.exists(parent)) {
            Files.createDirectories(parent);
          }
          // write file content
          Files.copy(zis, newFile, StandardCopyOption.REPLACE_EXISTING);
        }
        zipEntry = zis.getNextEntry();
      }
    }
  }


  /**
   * Protect against Zip Slip.
   */
  private static Path newFile(Path destDir, ZipEntry zipEntry) throws IOException {
    Path destFile = destDir.resolve(zipEntry.getName());
    Path normalizePath = destDir.normalize();
    if (!normalizePath.startsWith(destDir)) {
      throw new IOException("Bad zip entry: " + zipEntry.getName());
    }
    return destFile;
  }
}
