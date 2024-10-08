# Diplotype Comparison

This folder contains scripts for the following purposes:
1. For clinical calls, identify diplotypes that will always show up together as top-scored calls from PharmCAT.
2. For autogenerated tests, automatically check whether the autogenerated test results are plausible.

## Clinical calls

To obtain the pairs of diplotypes that always show up together as top-scored calls, run the following code:

```shell
python3 compare_diplotype_definition.py -c
```
This returns a tsv file named `predicted_pharmcat_calls.tsv` by default. It has only diplotypes that will show up together as possible calls for unphased data.

## Autogenerated tests

The same `compare_diplotype_definition.py` identifies all diplotypes that are defined by the same genotypes. These diplotypes may have different scores. In PharmCAT, diplotypes with a lower score will show up as alternative calls if the user chooses to output all possible calls. Nonetheless, a computational approach is needed to pas the autogenerated tests.

### 1. without missing positions
Run the following script to generate a reference file that contains all diplotypes defined by the same genotypes.

```shell
python3 compare_diplotype_definition.py
```
The default output is a file also named `predicted_pharmcat_calls.tsv`. The file lists
1. `expected`: expected call for unphased data
2. `actual`: the top-scored call that will be returned by PharmCAT
3. `alternative`: other lower-scored calls defined by the same genotypes as the `actual` call
4. `missing_positions`: a list of positions that are presumed missing
5. `gene`: gene symbols


### 2. with missing positions
Find possible alternative calls for every possible combination of missing positions. Note, to spped up the tool, the script only focuses on positions whose missing will lead to ambiguous calls.
```shell
python3 compare_diplotype_definition.py -m
```

Use the following command if you want to find alternative calls when specific positions are missing. At the moment, the script takes a VCF file and regards all positions in the VCF missing.
```shell
python3 compare_diplotype_definition.py -m -M missing_position.vcf
```

### 3. for phased data

Find alternative calls for phased data. This is specifically designed for phased data that have missing positions. Use `-M <missing_vcf>` if you want to specify the missing positions.

```shell
python3 compare_diplotype_definition.py -p
```



### 3. filter failed tests

Run the following code to pass autogenerated tests. Note that you need to update `test_file_pattern` to denote what autogenerated tests you are evaluating.
```shell
python3 filter_tests.py -f ~/Downloads/autogeneratedTestResults/autogenerated_test_report.tsv
```
