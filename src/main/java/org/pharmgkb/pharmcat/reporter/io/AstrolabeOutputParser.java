package org.pharmgkb.pharmcat.reporter.io;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import org.pharmgkb.pharmcat.reporter.model.AstrolabeCall;


/**
 * Class to parse an Astrolabe output file into {@link AstrolabeCall} objects
 *
 * @author Ryan Whaley
 */
public class AstrolabeOutputParser {

  private static final List<String> sf_geneWhitelist = ImmutableList.of("CYP2D6");

  @Nonnull
  public static List<AstrolabeCall> parse(@Nonnull Path astrolabePath) throws IOException {
    Preconditions.checkNotNull(astrolabePath);

    return Files.lines(astrolabePath)
        .filter(l -> !l.startsWith("#"))
        .map(AstrolabeCall::new)
        .filter(c -> sf_geneWhitelist.contains(c.getGene()))
        .collect(Collectors.toList());
  }
}
