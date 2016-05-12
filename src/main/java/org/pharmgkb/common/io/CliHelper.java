package org.pharmgkb.common.io;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import javax.annotation.Nonnull;
import com.google.common.collect.Lists;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;


/**
 * @author Mark Woon
 */
public class CliHelper {
  private static final String sf_verboseFlag = "verbose";
  private static final String sf_helpFlag = "help";
  private String m_name;
  private Options m_helpOptions = new Options();
  private Options m_options = new Options();
  private CommandLine m_commandLine;
  private String m_error;


  /**
   * Standard constructor.
   *
   * @param cls the class with the main() method
   */
  public CliHelper(Class cls) {

    m_name = cls.getSimpleName();

    Option opt = new Option("h", sf_helpFlag, false, "print this message");
    m_helpOptions.addOption(opt);
    m_options.addOption(opt);

    opt = new Option("v", sf_verboseFlag, false, "enable verbose output");
    m_helpOptions.addOption(opt);
    m_options.addOption(opt);
  }


  public CliHelper addOption(@Nonnull Option option) {

    if (option.getArgName().equals("h") || option.getArgName().equals("v")) {
      throw new IllegalArgumentException("-h and -v are reserved arguments");
    }
    m_helpOptions.addOption(option);
    m_options.addOption(option);
    return this;
  }

  /**
   * Add a boolean option (aka a flag).
   */
  public CliHelper addOption(@Nonnull String shortName, String longName, @Nonnull String description) {

    if (shortName.equals("h") || shortName.equals("v")) {
      throw new IllegalArgumentException("-h and -v are reserved arguments");
    }
    Option opt = Option.builder(shortName)
        .longOpt(longName)
        .desc(description)
        .hasArg(false)
        .build();
    m_helpOptions.addOption(opt);
    m_options.addOption(opt);
    return this;
  }

  /**
   * Add a required option that takes an argument.
   */
  public CliHelper addOption(@Nonnull String shortName, String longName, @Nonnull String description,
      boolean isOptionRequired, @Nonnull String argName) {
    return addOption(shortName, longName, description, isOptionRequired, argName, 1, true);
  }

  /**
   * Add a required option that takes arguments.
   *
   * @param numArgs 0 if argument(s) are optional, otherwise the number of expected arguments
   */
  public CliHelper addOption(@Nonnull String shortName, String longName, @Nonnull String description,
      boolean isOptionRequired, @Nonnull String argName, int numArgs, boolean argsAreRequired) {

    if (shortName.equals("h") || shortName.equals("v")) {
      throw new IllegalArgumentException("-h and -v are reserved arguments");
    }

    m_helpOptions.addOption(buildOption(shortName, longName, description, false, argName, numArgs, argsAreRequired));
    m_options.addOption(buildOption(shortName, longName, description, isOptionRequired, argName, numArgs, argsAreRequired));
    return this;
  }


  private Option buildOption(@Nonnull String shortName, String longName, @Nonnull String description,
      boolean isOptionRequired, @Nonnull String argName, int numArgs, boolean argsAreRequired) {

    Option.Builder optBuilder = Option.builder(shortName)
        .longOpt(longName)
        .desc(description)
        .argName(argName)
        .numberOfArgs(numArgs);

    if (argsAreRequired) {
      optBuilder.hasArg();
    } else {
      optBuilder.optionalArg(true);
    }
    // add non-require variant to help options
    if (isOptionRequired) {
      optBuilder.required();
    }
    return optBuilder.build();
  }


  /**
   * Parses arguments.
   *
   * @return true if parse completed and processing should continue, false if there are missing arguments or help was
   * requested
   */
  public boolean parse(String[] args) {

    try {
      CommandLineParser parser = new DefaultParser();
      // check for -h
      m_commandLine = parser.parse(m_helpOptions, args);
      if (isHelpRequested()) {
        printHelp();
        return false;
      }
      parser = new DefaultParser();
      m_commandLine = parser.parse(m_options, args);
      return true;

    } catch (org.apache.commons.cli.ParseException ex) {
      m_error = ex.getMessage();
      System.err.println(m_error);
      System.err.println();
      printHelp();
      return false;
    }
  }


  /**
   * Checks whether the specified option exists.
   */
  public boolean hasOption(String opt) {
    return m_commandLine.hasOption(opt);
  }

  /**
   * Gets the String value for the given option.
   */
  public String getValue(String opt) {
    return m_commandLine.getOptionValue(opt);
  }

  /**
   * Gets the String values for the given option.
   */
  public List<String> getValues(String opt) {
    return Lists.newArrayList(m_commandLine.getOptionValues(opt));
  }

  /**
   * Gets the int value for the given option.
   */
  public int getIntValue(String opt) {
    return Integer.parseInt(m_commandLine.getOptionValue(opt));
  }


  /**
   * Gets the value for the given option as a {@link File}.
   *
   * @param createIfNotExist if true and directory doesn't exist, create the directory;
   * otherwise, if false and directory doesn't exist, throw IllegalStateException
   * @return the directory
   * @throws IllegalStateException if directory doesn't exist
   */
  public @Nonnull Path getValidDirectory(String opt, boolean createIfNotExist) throws IOException {

    Path dir = Paths.get(getValue(opt));
    if (Files.exists(dir)) {
      if (Files.isDirectory(dir)) {
        return dir;
      }
      throw new IOException("Not a valid directory: " + dir);
    } else if (createIfNotExist) {
      Files.createDirectories(dir);
      return dir;
    }
    throw new IllegalStateException("No such directory: " + dir);
  }


  /**
   * Gets the value for the given option as a {@link Path}.
   *
   * @throws IllegalStateException if option was not specified
   */
  public @Nonnull Path getPath(@Nonnull String opt) {

    if (!hasOption(opt)) {
      throw new IllegalStateException("Missing option '" + opt + "'");
    }
    return Paths.get(getValue(opt));
  }

  /**
   * Gets the value for the given option as a {@link Path}, that must point to an existing file.
   *
   * @throws IllegalStateException if file doesn't exist
   */
  public @Nonnull Path getValidFile(@Nonnull String opt, boolean mustExist) {
    Path p = getPath(opt);
    if (!Files.exists(p)) {
      if (mustExist) {
        throw new IllegalStateException("File '" + p.toString() + "' does not exist");
      }
    } else {
      if (!Files.isRegularFile(p)) {
        throw new IllegalStateException("Not a file: '" + p.toString());
      }
    }
    return p;
  }


  /**
   * Gets remaining parameters.
   */
  public List getArguments() {
    return m_commandLine.getArgList();
  }


  /**
   * Gets whether to operate in verbose mode.
   */
  public boolean isVerbose() {
    return m_commandLine.hasOption(sf_verboseFlag);
  }


  /**
   * Checks whether the arguments were parsed successfully.
   */
  boolean hasError() {
    return m_error != null;
  }

  /**
   * Gets the error that occured while parsing arguments.
   */
  String getError() {
    return m_error;
  }


  /**
   * Gets whether help on command line arguments has been requested.
   */
  public boolean isHelpRequested() {
    return m_commandLine != null && m_commandLine.hasOption(sf_helpFlag);
  }


  /**
   * Prints the help message.
   */
  public void printHelp() {

    HelpFormatter formatter = new HelpFormatter();
    formatter.printHelp(m_name, m_options);
  }
}
