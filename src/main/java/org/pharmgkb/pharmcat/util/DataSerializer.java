package org.pharmgkb.pharmcat.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.concurrent.ThreadSafe;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.Sets;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;


/**
 * Serializer/Deserializer for data files.
 *
 * @author Mark Woon
 */
@ThreadSafe
public class DataSerializer {
  private static final Gson sf_gson = new GsonBuilder()
      .serializeNulls()
      .disableHtmlEscaping()
      .excludeFieldsWithoutExposeAnnotation()
      .setPrettyPrinting().create();
  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults().omitEmptyStrings();



  public void serializeToJson(@Nonnull Object data, @Nonnull Path jsonFile) throws IOException {
    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"), "Invalid format: %s does not end with .json", jsonFile);

    try (BufferedWriter writer = Files.newBufferedWriter(jsonFile, StandardCharsets.UTF_8)) {
      sf_gson.toJson(data, writer);
    }
  }


  public DefinitionFile deserializeDefinitionsFromJson(@Nonnull Path jsonFile) throws IOException {
    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"), "Invalid format: file name does not end with .json");
    Preconditions.checkArgument(Files.isRegularFile(jsonFile), "%s is not a file", jsonFile);

    try (BufferedReader reader = Files.newBufferedReader(jsonFile, StandardCharsets.UTF_8)) {
      DefinitionFile definitionFile = sf_gson.fromJson(reader, DefinitionFile.class);
      for (NamedAllele namedAllele : definitionFile.getNamedAlleles()) {
        namedAllele.initialize(definitionFile.getVariants());
      }
      return definitionFile;
    }
  }


  public Set<DefinitionExemption> deserializeExemptionsFromJson(@Nonnull Path jsonFile) throws IOException {
    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"), "Invalid format: %s does not end with .json", jsonFile);
    Preconditions.checkArgument(Files.isRegularFile(jsonFile), "%s is not a file", jsonFile);

    try (BufferedReader reader = Files.newBufferedReader(jsonFile, StandardCharsets.UTF_8)) {
      DefinitionExemption[] exemptions = sf_gson.fromJson(reader, DefinitionExemption[].class);
      return Sets.newHashSet(exemptions);
    }
  }

  public Set<DefinitionExemption> deserializeExemptionsFromTsv(@Nonnull Path tsvFile) throws IOException {
    Preconditions.checkNotNull(tsvFile);
    Preconditions.checkArgument(tsvFile.toString().endsWith(".tsv"), "Invalid format: %s does not end with .tsv", tsvFile);
    Preconditions.checkArgument(Files.isRegularFile(tsvFile), "%s is not a file", tsvFile);

    return Files.lines(tsvFile)
        .skip(1) // skip the header
        .map(line -> {
          String[] data = line.split("\t");
          String gene = data[0];
          SortedSet<String> ignoreAlleles = Sets.newTreeSet(sf_commaSplitter.splitToList(data[1]));
          boolean allHits = data.length > 2 && Boolean.parseBoolean(data[2]);
          boolean assumeReference = true;
          if (data.length > 3 && StringUtils.stripToEmpty(data[3]).equalsIgnoreCase("false")) {
            assumeReference = false;
          }
          return new DefinitionExemption(gene, ignoreAlleles, allHits, assumeReference);
        })
        .collect(Collectors.toSet());
  }


  public List<MessageAnnotation> deserializeMessagesFromTsv(@Nonnull Path tsvFile) throws IOException {
    Preconditions.checkNotNull(tsvFile);
    Preconditions.checkArgument(tsvFile.toString().endsWith(".tsv"), "Invalid format: %s does not end with .tsv", tsvFile);
    Preconditions.checkArgument(Files.isRegularFile(tsvFile), "%s is not a file", tsvFile);

    return Files.lines(tsvFile)
        .skip(1) // skip the header
        .map(MessageAnnotation::new)
        .collect(Collectors.toList());
  }
}
