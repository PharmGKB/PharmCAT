package org.pharmgkb.pharmcat.definition;

import java.io.IOException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * This class will take curator generated allele definition files and generate the "final" versions for consumption.
 *
 * @author Ryan Whaley
 * @author Alex Frase
 */
public class DefinitionTransformer {
  private static final Logger sf_logger = LoggerFactory.getLogger(DefinitionTransformer.class);


  public static void main(String[] args) {

    if (args.length != 2) {
      System.err.println("Missing input and output directory.");
      System.exit(1);
    }

    DefinitionTransformer formatter = new DefinitionTransformer();
    try {
      Path inDir = Paths.get(args[0]);
      if (!Files.exists(inDir)) {
        System.err.println("Input directory '" + args[0] + "' does not exist");
        System.exit(1);
      }
      if (!Files.isDirectory(inDir)) {
        System.err.println("Input directory '" + args[0] + "' is not a directory");
        System.exit(1);
      }

      Path outDir = Paths.get(args[0]);
      if (Files.exists(inDir) && !Files.isDirectory(inDir)) {
        System.err.println("Output directory '" + args[0] + "' is not a directory");
        System.exit(1);
      }
      Files.createDirectories(outDir);

      formatter.run(inDir, outDir);

    } catch (IOException ex) {
      ex.printStackTrace();
      System.exit(1);
    }
    System.exit(0);
  }


  /**
   * Does the work for stepping through the files and applying the format.
   * @throws IOException can occur due to Filesystem errors
   */
  private void run(Path inDir, Path outDir) throws IOException {

    try (DirectoryStream<Path> files = Files.newDirectoryStream(inDir, f -> f.toString().endsWith(".tsv"))) {
      for (Path file : files) {
        sf_logger.info("Parsing {}", file);
        CuratedDefinitionParser parser = new CuratedDefinitionParser(file);

        DefinitionFile definitionFile = parser.parse();
        if (parser.getWarnings().isEmpty()) {
          System.out.println("Warnings for " + file);
          parser.getWarnings().stream()
              .forEach(System.out::println);
        }

        Path outputFile = outDir.resolve(definitionFile.getGeneSymbol() + ".tsv");
        GeneratedDefinitionSerializer serializer = new GeneratedDefinitionSerializer();
        serializer.serializeToTsv(definitionFile, outputFile);
        sf_logger.info("Translated definition for {} to {}", definitionFile.getGeneSymbol(), outputFile);
      }
    }
  }
}
