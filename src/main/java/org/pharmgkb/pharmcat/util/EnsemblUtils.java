package org.pharmgkb.pharmcat.util;

import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.http.HttpEntity;
import org.apache.http.HttpHeaders;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClients;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;


/**
 * Utility class to download variant information from Ensembl
 *
 * @author Ryan Whaley
 */
public class EnsemblUtils {
  private static final String sf_variantUrl = "http://rest.ensembl.org/variation/human/%s";
  private static final String sf_targetAssembly = "GRCh38";
  private static final Gson sf_gson = new GsonBuilder()
      .serializeNulls()
      .disableHtmlEscaping()
      .excludeFieldsWithoutExposeAnnotation()
      .setPrettyPrinting().create();

  /**
   * Download variant data from Ensembl based on RSID and return the result as a populated VariantLocus.
   *
   * @param rsid the dbSNP RSID for the variant
   * @return a new VariantLocus with the data
   */
  @Nullable
  public static VariantLocus download(@Nonnull String rsid) {
    VariantLocus variantLocus = null;
    CloseableHttpClient httpclient = HttpClients.createDefault();
    HttpGet httpGet = new HttpGet(String.format(sf_variantUrl, rsid));
    httpGet.setHeader(HttpHeaders.ACCEPT, "application/json");
    try (CloseableHttpResponse response = httpclient.execute(httpGet)) {
      HttpEntity entity = response.getEntity();
      variantLocus = parse(entity.getContent());
    } catch (Exception e) {
      e.printStackTrace();
    }
    return variantLocus;
  }

  /**
   * Parse the given stream of JSON data into a VariantLocus with as much data populated as possible.
   *
   * Internally, this uses a {@link EnsemblResponse} class that's defined in this class
   * @param inputStream a stream of JSON text
   * @return a partially populated {@link VariantLocus} object
   */
  protected static VariantLocus parse(InputStream inputStream) {
    EnsemblResponse response = sf_gson.fromJson(new InputStreamReader(inputStream), EnsemblResponse.class);
    EnsemblMapping mapping = response.lookupMapping();
    VariantLocus variantLocus = new VariantLocus(mapping.getChrName(), mapping.getStart(), "g."+mapping.getStart());
    variantLocus.setRsid(response.getName());

    return variantLocus;
  }

  private class EnsemblMapping {
    private static final String sf_chrPrefix = "chr";

    @Expose
    @SerializedName("start")
    private int m_start;
    @Expose
    @SerializedName("seq_region_name")
    private String m_chrName;
    @Expose
    @SerializedName("assembly_name")
    private String m_assemlby;
    @Expose
    @SerializedName("allele_string")
    private String m_alleles;

    public int getStart() {
      return m_start;
    }

    public String getChrName() {
      return sf_chrPrefix + m_chrName;
    }

    public String getAssemlby() {
      return m_assemlby;
    }

    public String getAlleles() {
      return m_alleles;
    }
  }

  private class EnsemblResponse {

    @Expose
    @SerializedName("name")
    private String m_name;
    @Expose
    @SerializedName("mappings")
    private List<EnsemblMapping> m_mappings;
    @Expose
    @SerializedName("minor_allele")
    private String m_minorAllele;

    public String getName() {
      return m_name;
    }

    public List<EnsemblMapping> getMappings() {
      return m_mappings;
    }

    public String getMinorAllele() {
      return m_minorAllele;
    }

    public EnsemblMapping lookupMapping() {
      if (m_mappings == null || m_mappings.size() == 0) {
        return null;
      }

      return m_mappings.stream()
          .filter(m -> m.getAssemlby().equals(sf_targetAssembly))
          .findFirst().orElseThrow(RuntimeException::new);
    }

    public String lookupRefAllele() {
      EnsemblMapping mapping = lookupMapping();
      if (mapping == null) {
        return null;
      }

      String[] alleles = mapping.getAlleles().split("/");
      return Arrays.stream(alleles)
          .filter(a -> !a.equals(getMinorAllele()))
          .findFirst()
          .orElseThrow(RuntimeException::new);
    }
  }
}
