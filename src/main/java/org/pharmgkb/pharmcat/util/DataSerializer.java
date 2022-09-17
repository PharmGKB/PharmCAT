package org.pharmgkb.pharmcat.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.Sets;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;


/**
 * Serializer/Deserializer for data files.
 *
 * @author Mark Woon
 */
public class DataSerializer {
  public static final Gson GSON = new GsonBuilder()
      .serializeNulls()
      .disableHtmlEscaping()
      .excludeFieldsWithoutExposeAnnotation()
      .setPrettyPrinting()
      .create();
  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults().omitEmptyStrings();



  public void serializeToJson(Object data, Path jsonFile) throws IOException {
    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"), "Invalid format: %s does not end with .json", jsonFile);

    try (BufferedWriter writer = Files.newBufferedWriter(jsonFile, StandardCharsets.UTF_8)) {
      GSON.toJson(data, writer);
    }
  }


  public DefinitionFile deserializeDefinitionsFromJson(Path jsonFile) throws IOException {
    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"), "Invalid format: file name does not end with .json");
    Preconditions.checkArgument(Files.isRegularFile(jsonFile), "%s is not a file", jsonFile);

    try (BufferedReader reader = Files.newBufferedReader(jsonFile, StandardCharsets.UTF_8)) {
      DefinitionFile definitionFile = GSON.fromJson(reader, DefinitionFile.class);
      for (NamedAllele namedAllele : definitionFile.getNamedAlleles()) {
        namedAllele.initialize(definitionFile.getVariants());
      }
      return definitionFile;
    }
  }


  public Set<DefinitionExemption> deserializeExemptionsFromJson(Path jsonFile) throws IOException {
    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"), "Invalid format: %s does not end with .json", jsonFile);
    Preconditions.checkArgument(Files.isRegularFile(jsonFile), "%s is not a file", jsonFile);

    try (BufferedReader reader = Files.newBufferedReader(jsonFile, StandardCharsets.UTF_8)) {
      DefinitionExemption[] exemptions = GSON.fromJson(reader, DefinitionExemption[].class);
      return Sets.newHashSet(exemptions);
    }
  }

  public Set<DefinitionExemption> deserializeExemptionsFromTsv(Path tsvFile) throws IOException {
    Preconditions.checkNotNull(tsvFile);
    Preconditions.checkArgument(tsvFile.toString().endsWith(".tsv"), "Invalid format: %s does not end with .tsv", tsvFile);
    Preconditions.checkArgument(Files.isRegularFile(tsvFile), "%s is not a file", tsvFile);

    SortedSet<DefinitionExemption> exemptions = new TreeSet<>();
    try (VcfHelper vcfHelper = new VcfHelper();
         BufferedReader reader = Files.newBufferedReader(tsvFile)) {
      // skip the header
      reader.readLine();
      String line;
      while ((line  = reader.readLine()) != null) {
        String[] data = line.split("\t");
        String gene = data[0];
        final SortedSet<VariantLocus> ignoreLoci = new TreeSet<>();
        if (data.length > 1) {
          for (String rsid : sf_commaSplitter.splitToList(data[1])) {
            ignoreLoci.add(vcfHelper.rsidToVariantLocus(rsid));
          }
        }
        final SortedSet<VariantLocus> extraLoci = new TreeSet<>();
        if (data.length > 2) {
          for (String rsid : sf_commaSplitter.splitToList(data[2])) {
            extraLoci.add(vcfHelper.rsidToVariantLocus(rsid));
          }
        }
        SortedSet<String> ignoreAlleles = null;
        if (data.length > 3) {
          ignoreAlleles = Sets.newTreeSet(sf_commaSplitter.splitToList(data[3]));
        }
        Boolean allHits = null;
        if (data.length > 4 && StringUtils.stripToNull(data[4]) != null) {
          allHits = Boolean.parseBoolean(data[4]);
        }
        exemptions.add(new DefinitionExemption(gene, ignoreLoci, extraLoci, ignoreAlleles, allHits));
      }
    }
    return exemptions;
  }


  public List<MessageAnnotation> deserializeMessagesFromTsv(Path tsvFile) throws IOException {
    Preconditions.checkNotNull(tsvFile);
    Preconditions.checkArgument(tsvFile.toString().endsWith(".tsv"), "Invalid format: %s does not end with .tsv", tsvFile);
    Preconditions.checkArgument(Files.isRegularFile(tsvFile), "%s is not a file", tsvFile);

    List<MessageAnnotation> messageAnnotations = new ArrayList<>();
    try (BufferedReader reader = Files.newBufferedReader(tsvFile)) {
      // skip the header
      reader.readLine();
      String line;
      int x = 1;
      while ((line = reader.readLine()) != null) {
        x += 1;
        if (StringUtils.isBlank(line)) {
          continue;
        }
        MessageAnnotation ma = new MessageAnnotation(line);
        if (ma.getName() == null) {
          throw new IllegalStateException("Row " + x + ": Missing name");
        }
        messageAnnotations.add(ma);
      }
    }
    return messageAnnotations;
  }
}
