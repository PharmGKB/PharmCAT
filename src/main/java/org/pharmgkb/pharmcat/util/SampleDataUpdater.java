package org.pharmgkb.pharmcat.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveOutputStream;
import org.apache.commons.io.IOUtils;
import org.jspecify.annotations.Nullable;


/**
 * This tool updates the sample data in the {@code docs/examples} directory.
 *
 * @author Mark Woon
 */
public class SampleDataUpdater {

  public static void main(String[] args) {

    Path pharmcatDir = null;
    if (args.length != 1) {
      if (System.getenv("PHARMCAT") != null) {
        pharmcatDir = Paths.get(System.getenv("PHARMCAT"));
      } else {
        System.out.println("Expecting a single argument - the PharmCAT directory");
        System.exit(1);
      }
    } else {
      pharmcatDir = Paths.get(args[0]);
    }
    if (!Files.isDirectory(pharmcatDir)) {
      System.out.println("Not a valid directory: " + pharmcatDir);
      System.exit(1);
    }

    Path positionsFile = pharmcatDir.resolve("pharmcat_positions.vcf");
    if (!Files.exists(positionsFile)) {
      System.out.println("Example file does not exist: " + positionsFile);
      System.exit(1);
    }

    Path examplesDir = pharmcatDir.resolve("docs/examples");
    if (!Files.exists(examplesDir)) {
      System.out.println("Examples directory does not exist: " + examplesDir);
      System.exit(1);
    }

    try {
      Path sample1File = examplesDir.resolve("pharmcat.example.vcf");
      try (BufferedReader reader = Files.newBufferedReader(positionsFile);
           PrintWriter writer = new PrintWriter(Files.newBufferedWriter(sample1File))) {
        String line;
        while ((line = reader.readLine()) != null) {
          line = replace(line, null, "FORMAT\tPharmCAT", "FORMAT\tSample_1");
          writer.println(line);
        }
      }

      Path sample2File = examplesDir.resolve("pharmcat.example2.vcf");
      try (BufferedReader reader = Files.newBufferedReader(positionsFile);
           PrintWriter writer = new PrintWriter(Files.newBufferedWriter(sample2File))) {
        String line;
        while ((line = reader.readLine()) != null) {
          line = replace(line, new String[][] {
              new String[] {null, "FORMAT\tPharmCAT", "FORMAT\tSample_2"},
              new String[] {"rs12769205", "0/0", "1/1"},
              new String[] {"rs4244285", "0/0", "1/1"},
              new String[] {"rs3758581", "0/0", "1/1"},
              new String[] {"rs3745274", "0/0", "0/1"},
              new String[] {"rs2279343", "0/0", "0/1"}
          });
          writer.println(line);
        }
      }

      Path multisampleFile = examplesDir.resolve("multisample.vcf");
      try (BufferedReader reader1 = Files.newBufferedReader(sample1File);
           BufferedReader reader2 = Files.newBufferedReader(sample2File);
           PrintWriter writer = new PrintWriter(Files.newBufferedWriter(multisampleFile))) {
        String line1;
        while ((line1 = reader1.readLine()) != null) {
          String line2 = reader2.readLine();
          if (line1.startsWith("#")) {
            if (line1.startsWith("##bcftools")) {
              continue;
            } else if (line1.startsWith("##source=")) {
              line1 = "##source=PharmCAT examples";
            } else if (line1.endsWith("FORMAT\tSample_1")) {
              line1 += "\tSample_2";
            }
            writer.println(line1);
            continue;
          }

          String[] data = line2.split("\t");
          writer.println(line1 + "\t" + data[data.length - 1]);
        }
      }


      Path splitListFile = examplesDir.resolve("split_vcf_list.txt");
      PrintWriter writer = null;
      List<Path> chrFiles = new ArrayList<>();
      try (BufferedReader reader1 = Files.newBufferedReader(multisampleFile);
           PrintWriter listWriter = new PrintWriter(Files.newBufferedWriter(splitListFile))) {
        String line;
        StringBuilder headers = new StringBuilder();
        String chr = null;
        while ((line = reader1.readLine()) != null) {
          if (line.startsWith("#")) {
            headers.append(line).append("\n");
            continue;
          }
          String[] data = line.split("\t");
          if (chr == null) {
            chr = data[0];
            writer = getNewChrWriter(examplesDir, listWriter, chr, chrFiles);
            writer.print(headers);
          } else if (!chr.equals(data[0])) {
            chr = data[0];
            writer.close();
            writer = getNewChrWriter(examplesDir, listWriter, chr, chrFiles);
            writer.print(headers);
          }
        }
        if (writer != null) {
          writer.close();
        }
      }

      Path splitFileTar = examplesDir.resolve("split_vcf.tar");
      try (OutputStream tarFileStream = Files.newOutputStream(splitFileTar);
           TarArchiveOutputStream tarStream = new TarArchiveOutputStream(tarFileStream)) {
        for (Path file : chrFiles) {
          addToArchive(tarStream, file);
        }
      }
      for (Path file : chrFiles) {
        Files.delete(file);
      }


      System.out.println("Done.");

    } catch (Exception ex) {
      //noinspection CallToPrintStackTrace
      ex.printStackTrace();
      System.exit(1);
    }
  }


  private static PrintWriter getNewChrWriter(Path dir, PrintWriter listWriter, String chr, List<Path> chrFiles)
      throws IOException {
    Path file = dir.resolve(chr + ".vcf");
    chrFiles.add(file);
    listWriter.println(chr + ".vcf");
    return new PrintWriter(Files.newBufferedWriter(file));
  }

  private static void addToArchive(TarArchiveOutputStream out, Path file) throws IOException {
    TarArchiveEntry entry = new TarArchiveEntry(file, file.getFileName().toString());
    // set time/user/group to 0 to make the tar file deterministic (based on content alone, so get the same md5 hash)
    entry.setModTime(0);
    entry.setUserId(0);
    entry.setGroupId(0);
    out.putArchiveEntry(entry);
    try (InputStream fis = Files.newInputStream(file)) {
      IOUtils.copy(fis, out);
    }
    out.closeArchiveEntry();
  }

  private static String replace(String line, String[][] modifications) {
    for (String[] mod : modifications) {
      String l = replace(line, mod[0], mod[1], mod[2]);
      if (!l.equals(line)) {
        return l;
      }
    }
    return line;
  }

  private static String replace(String line, @Nullable String marker, String from, String to) {
    if (marker == null) {
      marker = from;
    }
    if (line.contains(marker)) {
      line = line.replaceAll(from, to);
    }
    return line;
  }
}
