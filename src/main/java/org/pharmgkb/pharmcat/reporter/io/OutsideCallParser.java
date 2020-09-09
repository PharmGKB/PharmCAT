package org.pharmgkb.pharmcat.reporter.io;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;


/**
 * Class to parse outside-call text files into {@link OutsideCall} objects. This was created with Astrolabe in mind but 
 * can be used for any outside calls. This parser will ignore any lines starting with "#" or blank lines and interpret 
 * the rest as TSV formatted lines with the first field as gene symbol and the second field as diplotype call.
 *
 * @author Ryan Whaley
 */
public class OutsideCallParser {

  private static final List<String> sf_geneWhitelist = ImmutableList.of("CYP2D6");

  @Nonnull
  public static List<OutsideCall> parse(Path filePath) throws IOException {
    Preconditions.checkNotNull(filePath);

    return Files.lines(filePath)
        .filter(l -> !l.startsWith("#"))
        .map(OutsideCall::new)
        .filter(c -> sf_geneWhitelist.contains(c.getGene()))
        .collect(Collectors.toList());
  }
}
