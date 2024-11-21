# JSON-to-TSV Converter

This folder contains scripts that convert PharmCAT JSON files to a tabular TSV file for the ease of downstream analyses.

## How to run

Suppose you are running the `json2tsv_pharmcat.py` script from the folder where the script locates.
```console
# python3 json2tsv_pharmcat.py -i /path/to/pharmcat/json/files/
```
**Mandatory** argument: `-i`.

-i `</path/to/pharmcat/json/>`
: Path to the directory where PharmCAT Named Allele Matcher and Phenotyper JSON files are **BOTH** locates. By default, the script looks for files ending with `*.match.json` and `*.phenotype.json`.

#### Optional Arguments

-S `<txt_file>` <span class="altArg"><br />or --sample-file `<txt_file>`</span>
: A file of samples whose JSON files will be processed and included in the TSV. The file should contain one sample per line.

-m `<pattern>` <span class="altArg"><br />or --matcher-json-pattern `<pattern>`</span>
: Regular-expression-styled filename pattern for the PharmCAT Named Allele Matcher JSON files, by default `*.match.json`.

-p `<pattern>` <span class="altArg"><br />or --phenotyper-json-pattern `<pattern>`</span>
: Regular-expression-styled filename pattern for the PharmCAT Phenotyper JSON files, by default `*.phenotype.json`.

-a `</path/to/*_translation.json>` <span class="altArg"><br />or --allele-definition-files `</path/to/*_translation.json>`</span>
: Path to allele definition JSON files. The script will use the PharmCAT allele definition JSON files by default. If a path is provided, the script assumes the filename following the pattern `*_translation.json`. Include pattern matching strings if you are using custom allele definition JSON files.

-gs `CPIC` <span class="altArg"><br />or --guideline-source `CPIC/DPWG`</span>
: Guideline source to extract, default = `CPIC`. 

-g `<genes>` <span class="altArg"><br />or --genes `<genes>>`</span>
: List of genes to be processed. Separate by comma, `gene1,gene2,etc.`

-o `<dir>` <span class="altArg"><br />or --output-dir `<dir>`</span>
: Directory to save the TSV to. Default is the directory of the input.

-bf `<name>` <span class="altArg"><br />or --base-filename `<name>`</span>
: Prefix of the output TSV file (without file extension). Default is `pharmcat_json2tsv`.


-c <span class="altArg"><br />or --concurrent-mode</span>
: Enable concurrent mode. This defaults to using one CPU core if `-cp` is not specified. 

-cp `<num processes>` <span class="altArg"><br />or --max-concurrent-processes `<num processes>`</span>
: The maximum number of processes to use if concurrent mode is enabled.

-v <span class="altArg"><br />or --verbose</span>
: Print verbose messages.


## List of genes by guideline sources

The following is the default list of genes that the script extracts based on the guideline source `-gs`. You can specify the list of genes to extract by `-g`.

### CPIC

```console
ABCG2, CACNA1S, CFTR, CYP2B6, CYP2C9, 
CYP2C19, CYP3A5, CYP4F2, DPYD, G6PD, 
IFNL3, NAT2, NUDT15, RYR1, SLCO1B1, 
TPMT, UGT1A1, VKORC1, CYP2D6, HLA-A, 
HLA-B, MT-RNR1
```

### DPWG

```console
ABCG2, CYP2B6, CYP2C9, CYP2C19, CYP3A4, 
CYP3A5, DPYD, NUDT15, SLCO1B1, TPMT,
UGT1A1, VKORC1, CYP2D6
```

## Samples
The script identifies samples based on filenames of the input Matcher and Phenotyper JSON files. PharmCAT JSON files are named <file_basename>.<sample_name>.<match/phenotype>.json by default. The script assumes the sample names are located in the third last field of a filename and extracts the list of samples from Matcher and Phenotyper JSON filenames separately. 


Only samples that have both Matcher and Phenotyper JSON files will be passed to the downstream analysis. We need both Matcher and Phenotyper JSON files for comprehensive information about PGx calls and allele-defining genotypes. 

## JSON fields for PharmCAT data

* `componentHaplotypes` is a field in the Matcher JSON for combination alleles.

### All genes except _DPYD_ and _RYR1_
| TSV column                  | source                                   | JSON field, if application                                     |
|-----------------------------|------------------------------------------|----------------------------------------------------------------|
| sample                      | JSON filename                            | N/A                                                            |
| gene                        | Guideline source `-gs` or gene list `-g` | N/A                                                            |
| phenotype                   | Phenotyper JSON                          | phenotypes                                                     |
| activity_score              | Phenotyper JSON                          | activityScore, if applicable to the gene                       |
| diplotype                   | Phenotyper JSON                          | recommendationDiplotypes -> label                              |
| haplotype_1                 | Phenotyper JSON                          | recommendationDiplotypes -> allele1 -> name                    |
| haplotype_2                 | Phenotyper JSON                          | recommendationDiplotypes -> allele2 -> name (if available)     |
| haplotype_1_functions       | Phenotyper JSON                          | recommendationDiplotypes -> allele1 -> function                |
| haplotype_2_functions       | Phenotyper JSON                          | recommendationDiplotypes -> allele2 -> function (if available) |
| haplotype_1_variants        | Matcher JSON                             | diplotypes -> haplotype1 -> (componentHaplotypes) -> sequences |
| haplotype_2_variants        | Matcher JSON                             | diplotypes -> haplotype2 -> (componentHaplotypes) -> sequences |
| missing_positions           | Matcher JSON                             | missingPositions                                               |
| uncallable_haplotypes       | Matcher JSON                             | uncallableHaplotypes                                           |


### _DPYD_ and _RYR1_

#### Effectively phased
| TSV column                  | source           | JSON field, if application                                                |
|-----------------------------|------------------|---------------------------------------------------------------------------|
| phenotype                   | Phenotyper JSON  | recommendationDiplotypes -> phenotypes                                    |
| activity_score              | Phenotyper JSON  | recommendationDiplotypes -> activityScore                                 |
| diplotype                   | Phenotyper JSON  | sourceDiplotypes -> label                                                 |
| dpyd_ryr1_variants          | Phenotyper JSON  | matcherComponentHaplotypes -> allele1 -> name                             |
| dpyd_ryr1_variant_functions | Phenotyper JSON  | matcherComponentHaplotypes -> allele1 -> function                         |
| dpyd_ryr1_variant_genotypes | Matcher JSON     | diplotypes -> haplotype1/haplotype2 -> (componentHaplotypes) -> sequences |


#### Not effectively phased

| TSV column                  | source           | JSON field, if application                |
|-----------------------------|------------------|-------------------------------------------|
| phenotype                   | Phenotyper JSON  | recommendationDiplotypes -> phenotypes    |
| activity_score              | Phenotyper JSON  | recommendationDiplotypes -> activityScore |
| diplotype                   | Phenotyper JSON  | N/A                                       |
| dpyd_ryr1_variants          | Phenotyper JSON  | sourceDiplotypes -> allele1 -> name       |
| dpyd_ryr1_variant_functions | Phenotyper JSON  | sourceDiplotypes -> allele1 -> function   |
| dpyd_ryr1_variant_genotypes | Matcher JSON     | haplotypeMatches -> sequences             |

