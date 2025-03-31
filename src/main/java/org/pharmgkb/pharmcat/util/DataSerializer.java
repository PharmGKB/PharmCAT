package org.pharmgkb.pharmcat.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.reflect.Type;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.Sets;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.reflect.TypeToken;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
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
      .registerTypeAdapter(Date.class, new GsonDateAdapter())
      .setPrettyPrinting()
      .create();
  private static final Pattern sf_rsidPattern = Pattern.compile("rs\\d+");
  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults().omitEmptyStrings();
  private static final Splitter sf_semicolonSplitter = Splitter.on(";").trimResults().omitEmptyStrings();



  public static void serializeToJson(Object data, Path jsonFile) throws IOException {
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


  public SortedMap<String, DefinitionExemption> deserializeExemptionsFromJson(Path jsonFile) throws IOException {
    Preconditions.checkNotNull(jsonFile);
    Preconditions.checkArgument(jsonFile.toString().endsWith(".json"), "Invalid format: %s does not end with .json", jsonFile);
    Preconditions.checkArgument(Files.isRegularFile(jsonFile), "%s is not a file", jsonFile);

    try (BufferedReader reader = Files.newBufferedReader(jsonFile, StandardCharsets.UTF_8)) {
      Type mapType = new TypeToken<SortedMap<String, DefinitionExemption>>() {}.getType();
      return GSON.fromJson(reader, mapType);
    }
  }


  SortedMap<String, DefinitionExemption> deserializeExemptionsFromTsv(Path exemptionsTsvFile,
      @Nullable Path unphasedPrioritiesTsvFile) throws IOException {
    Preconditions.checkNotNull(exemptionsTsvFile);
    Preconditions.checkArgument(exemptionsTsvFile.toString().endsWith(".tsv"),
        "Invalid format: %s does not end with .tsv", exemptionsTsvFile);
    Preconditions.checkArgument(Files.isRegularFile(exemptionsTsvFile),
        "%s is not a file", exemptionsTsvFile);

    SortedMap<String, DefinitionExemption> exemptions = new TreeMap<>();
    try (VcfHelper vcfHelper = new VcfHelper();
         BufferedReader reader = Files.newBufferedReader(exemptionsTsvFile)) {
      // skip the header
      reader.readLine();
      String line;
      while ((line  = reader.readLine()) != null) {
        String[] data = line.split("\t");
        String gene = data[0];
        final SortedSet<VariantLocus> requiredLoci = new TreeSet<>();
        if (data.length > 1) {
          parseLoci(vcfHelper, data[1], requiredLoci);
        }
        final SortedSet<VariantLocus> ignoreLoci = new TreeSet<>();
        if (data.length > 2) {
          parseLoci(vcfHelper, data[2], ignoreLoci);
        }
        final SortedSet<VariantLocus> extraLoci = new TreeSet<>();
        if (data.length > 3) {
          parseLoci(vcfHelper, data[3], extraLoci);
        }
        SortedSet<String> ignoreAlleles = null;
        if (data.length > 4) {
          ignoreAlleles = Sets.newTreeSet(sf_commaSplitter.splitToList(data[4]));
        }
        List<String> amp1Alleles = null;
        if (data.length > 5) {
          amp1Alleles = sf_semicolonSplitter.splitToList(data[5]);
        }
        SortedSet<Long> requiredPositions = requiredLoci.stream()
            .map(VariantLocus::getPosition)
            .collect(Collectors.toCollection(TreeSet::new));
        exemptions.put(gene, new DefinitionExemption(gene, requiredPositions, ignoreLoci, extraLoci, ignoreAlleles,
            amp1Alleles));
      }
    }

    if (unphasedPrioritiesTsvFile != null) {
      try (BufferedReader reader = Files.newBufferedReader(unphasedPrioritiesTsvFile)) {
        // skip the header
        reader.readLine();
        String line;
        while ((line  = reader.readLine()) != null) {
          String[] data = line.split("\t");
          String gene = data[0];
          if (!gene.toUpperCase().equals(gene)) {
            throw new IllegalArgumentException("Gene name is not all uppercase: " + gene);
          }

          DefinitionExemption exemption = exemptions.computeIfAbsent(gene,
              g -> new DefinitionExemption(gene, null, null, null, null, null));

          String pick = data[1];
          SortedSet<String> dips = new TreeSet<>();
          //noinspection ManualArrayToCollectionCopy
          for (int x = 1; x < data.length; x++) {
            dips.add(data[x]);
          }
          exemption.addUnphasedDiplotypePriority(dips, pick);
        }
      }
    }
    return exemptions;
  }

  private void parseLoci(VcfHelper vcfHelper, String data, SortedSet<VariantLocus> loci) throws IOException {
    if (StringUtils.stripToNull(data) == null) {
      return;
    }
    for (String var : sf_commaSplitter.splitToList(data)) {
      if (sf_rsidPattern.matcher(var).matches()) {
        loci.add(vcfHelper.rsidToVariantLocus(var));
      } else {
        loci.add(vcfHelper.hgvsToVariantLocus(var));
      }
    }
  }

  public List<MessageAnnotation> deserializeMessagesFromTsv(Path tsvFile) throws IOException {
    Preconditions.checkNotNull(tsvFile);
    Preconditions.checkArgument(tsvFile.toString().endsWith(".tsv"), "Invalid format: %s does not end with .tsv", tsvFile);
    Preconditions.checkArgument(Files.isRegularFile(tsvFile), "%s is not a file", tsvFile);

    Set<String> keys = new HashSet<>();
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
        try {
          MessageAnnotation ma = new MessageAnnotation(line);
          if (ma.getName() == null) {
            throw new IllegalStateException("Row " + x + ": Missing name");
          }
          if (!keys.add(ma.getName())) {
            throw new IllegalStateException("Row " + x + ": Duplicate name '" + ma.getName() + "'");
          }
          messageAnnotations.add(ma);
        } catch (IllegalArgumentException ex) {
          throw new IllegalStateException("Row " + x + ": " + ex.getMessage(), ex);
        }
      }
    }
    String[] requiredKeys = new String[] {
        MessageHelper.MSG_COMBO_NAMING,
        MessageHelper.MSG_COMBO_UNPHASED,
        MessageHelper.MSG_CYP2D6_MODE,
        MessageHelper.MSG_CYP2D6_NOTE,
        MessageHelper.MSG_DPYD_HAPB3_EXONIC_ONLY,
        MessageHelper.MSG_DPYD_HAPB3_INTRONIC_MISMATCH_EXONIC,
        MessageHelper.MSG_MULTI_CALL,
        MessageHelper.MSG_MULTI_SCORE,
        MessageHelper.MSG_OUTSIDE_CALL,
    };
    for (String key : requiredKeys) {
      if (!keys.contains(key)) {
        throw new IllegalStateException("Missing static key: " + key);
      }
    }
    return messageAnnotations;
  }
}
