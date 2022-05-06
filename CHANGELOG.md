## [1.8.0](https://github.com/PharmGKB/PharmCAT/compare/v1.7.0...v1.8.0) (2022-05-06)


### Features

* update to latest PharmVar (5.1.14) data for CYP2C9 ([2d9de55](https://github.com/PharmGKB/PharmCAT/commit/2d9de55dd94274e6fd451e017526a33a5471843c))


### Bug Fixes

* fix VcfReaderTest for updated wording ([870b3b9](https://github.com/PharmGKB/PharmCAT/commit/870b3b9af9b8e95c515296cc72c49d1a7a311a59))
* **data:** update to CPIC 1.17.1 and PharmVar 5.1.14 ([d205472](https://github.com/PharmGKB/PharmCAT/commit/d20547215f0d13a221dfc393d4f122a1ba49a35a))
* **preprocessor:** clarify the handling of mismatched VCF entries ([b284637](https://github.com/PharmGKB/PharmCAT/commit/b284637d1928681aa993d2b02167b711a012d111))

## [1.7.0](https://github.com/PharmGKB/PharmCAT/compare/v1.6.0...v1.7.0) (2022-04-26)


### Features

* add a tracking event to the download button ([04f2062](https://github.com/PharmGKB/PharmCAT/commit/04f206211171eb68ef7e25efd10e5aff0d1da02c))
* add umami analytics to the pharmcat.org site ([82d8cd8](https://github.com/PharmGKB/PharmCAT/commit/82d8cd8215db45b21d4d38cab09b7958b365d324))
* CYP2D6 code and documentation for StellarPGx ([0b400ac](https://github.com/PharmGKB/PharmCAT/commit/0b400ac09832eabe8c229f6982542a57ee946665))
* update to v1.17 of CPIC and v5.1.12 of PharmVar ([0dcf2c9](https://github.com/PharmGKB/PharmCAT/commit/0dcf2c93d71cb5112fbbe20a46642da2024b7ff5))
* **preprocessor:** handle alt=<*> for INDELs ([776f32a](https://github.com/PharmGKB/PharmCAT/commit/776f32abea5bbb6489f552a48cae75a4743bd568))
* **reporter:** add a "test mode" to the Reporter ([7799d91](https://github.com/PharmGKB/PharmCAT/commit/7799d9118a6ea11dad398353382226525a287e01))
* **reporter:** add outside phenotype calls and support allele status genes ([d528d7a](https://github.com/PharmGKB/PharmCAT/commit/d528d7a92343ea20a8808aef4cd64729b97dac87))


### Bug Fixes

* fix command line arg for supplying your own named alleles ([4a8c84f](https://github.com/PharmGKB/PharmCAT/commit/4a8c84fa077d869b619e8e603509a7f2dd73fe6a))
* update site config ([07d23cc](https://github.com/PharmGKB/PharmCAT/commit/07d23cc0dea146a22db26d79f6f73ba52b185758))
* **namedallelematcher:** improve warnings about structural variations ([63903df](https://github.com/PharmGKB/PharmCAT/commit/63903df00f8c1fc744920a0571d3e35ad416f0bc))
* **preprocessor:** handle alt=<*> (unspecific alleles) for SNPs ([9c70001](https://github.com/PharmGKB/PharmCAT/commit/9c70001b6c8b393c46a977952fcf53358e92dd3e))
* **reporter:** fix handling of G6PD ([cbcfb05](https://github.com/PharmGKB/PharmCAT/commit/cbcfb050ac15269a319af1cb23025b278c7cbbdf))
* **reporter:** fix test outside call files for new syntax ([39ed203](https://github.com/PharmGKB/PharmCAT/commit/39ed20328a2ab6efd083aa5ee8cd80eddc5f251a))
* **reporter:** hide genes that are not reportable in the genotype table of the report ([51a1bf7](https://github.com/PharmGKB/PharmCAT/commit/51a1bf7ffc4e51451ebe7b3909a2e45fbaf6e82f))

## [1.6.0](https://github.com/PharmGKB/PharmCAT/compare/v1.5.1...v1.6.0) (2022-03-31)


### Features

* pharmcat.org website redesign and new content ([48cd154](https://github.com/PharmGKB/PharmCAT/commit/48cd15486bd342c96473e633f19ce1a735eb647b))
* **website:** tweak the website styling and layout ([b52b7d6](https://github.com/PharmGKB/PharmCAT/commit/b52b7d6f2bdd820cf506aa774c4ccd7beb5c1da3))


### Bug Fixes

* add more clarifying info to TestVcfBuilder alt allele mismatch message ([63d0dfd](https://github.com/PharmGKB/PharmCAT/commit/63d0dfd707d52317540d1a4b9095a8d190f3b69f))
* fix bundle call for running local jekyll ([a1034a7](https://github.com/PharmGKB/PharmCAT/commit/a1034a7c61494e2ba946775e3f0fdce37e0da1b2))
* fix CLI arg description for Phenotyper ([f31085c](https://github.com/PharmGKB/PharmCAT/commit/f31085c279c4ed163c665b6af0ebee438076bbc7))
* fix jekyll front matter on disclaimers page ([a2ecb79](https://github.com/PharmGKB/PharmCAT/commit/a2ecb790559b871d60767e15e9e40cd64f41f525))
* fix permissions on VCF scripts ([3d5cb34](https://github.com/PharmGKB/PharmCAT/commit/3d5cb3403f3755096bd4d0020d24077233093b65)), closes [#89](https://github.com/PharmGKB/PharmCAT/issues/89)
* **data:** update to CPIC 1.16 ([5ee8d37](https://github.com/PharmGKB/PharmCAT/commit/5ee8d37e7b4560c6db1c74eb97f2855634c553ce))
* **datamanager:** correctly loading definitions in dataManager when skipping alleles ([2687635](https://github.com/PharmGKB/PharmCAT/commit/26876359a58944a299be618e4e2760c85e3f6b53))
* **datamanager:** improve error when encounter external service error ([93de324](https://github.com/PharmGKB/PharmCAT/commit/93de324bcbccf5f45ebb7eabf2afdef2bee91a38))
* **docker:** make sure scripts are executable ([dfc080f](https://github.com/PharmGKB/PharmCAT/commit/dfc080f305da4a8e7e26ff78b38ee180ccefd3c0)), closes [#89](https://github.com/PharmGKB/PharmCAT/issues/89)
* **docker:** use docker 3.9 to work around numpy incompatibility ([44cf6d1](https://github.com/PharmGKB/PharmCAT/commit/44cf6d16552132c156bf44f4aed72dc3778af3fc))
* **namedallelematcher:** discard position with AD field ([6b717db](https://github.com/PharmGKB/PharmCAT/commit/6b717db24f9e64ab06f943958665dab7c3b71be7))
* **namedallelematcher:** do not discard position if novel ALT allele is found ([ba48472](https://github.com/PharmGKB/PharmCAT/commit/ba48472e4f3d44281f6efecaf0587f7acefe7972))
* **namedallelematcher:** handle `.` AD value ([2591ae9](https://github.com/PharmGKB/PharmCAT/commit/2591ae94fdcda24fdd88ec9bd9a4061c0def46a9))
* **namedallelematcher:** improve AD field handling (only catch reference overlap) ([cb152c6](https://github.com/PharmGKB/PharmCAT/commit/cb152c64399299bad9ecf2edd98abbefef2fe559))
* **namedallelematcher:** improve checking for unexpected ALT ([f74b309](https://github.com/PharmGKB/PharmCAT/commit/f74b309202affb07f0c0d7c35ff18f5a4b801346))
* **preprocessor:** improved way to identify block gVCF ([b57a6e2](https://github.com/PharmGKB/PharmCAT/commit/b57a6e2b1bb31cbaa42e5e08a5e883e83d35c0cc)), closes [#79](https://github.com/PharmGKB/PharmCAT/issues/79)
* **preprocessor:** sort VCF in a way that the non-PGx variants that occur at the PGx positions will appear after the line of the PGx variants ([3a2be37](https://github.com/PharmGKB/PharmCAT/commit/3a2be374a4973708aaa50dafae95371e7b0c1400)), closes [#95](https://github.com/PharmGKB/PharmCAT/issues/95)
* **reporter:** fix total gene count in report genotypes table ([02fd57d](https://github.com/PharmGKB/PharmCAT/commit/02fd57d4d3a30388c07847d22bc87669c0329d9a))
* **reporter:** make DrugLinks comparable to sort properly in final output ([501db22](https://github.com/PharmGKB/PharmCAT/commit/501db222b196b01100bb16c7cf8367d4c1ef0bab))
* **reporter:** update more PharmCATTest tests to new wrapper ([89e5354](https://github.com/PharmGKB/PharmCAT/commit/89e53548584f5ba8c7bedf2ec876954fdfb1aab0))


### Reverts

* update TPMT test for new calling rules ([63f4775](https://github.com/PharmGKB/PharmCAT/commit/63f4775e20cea117fa8ea8e61ea9dc3f31e8a644))

### [1.5.1](https://github.com/PharmGKB/PharmCAT/compare/v1.5.0...v1.5.1) (2022-02-22)


### Bug Fixes

* **data:** add ABCG2 ([92083db](https://github.com/PharmGKB/PharmCAT/commit/92083dbbe4bfdf20b32aa312960224b7de0cfda4))
* **reporter:** hide gene in the genotype summary table if no data present ([9db9de1](https://github.com/PharmGKB/PharmCAT/commit/9db9de16d72885b155616d85b09b731ac62f10bc))

## [1.5.0](https://github.com/PharmGKB/PharmCAT/compare/v1.4.0...v1.5.0) (2022-02-19)


### Features

* **reporter:** show Matcher variant warnings in the final report ([c247c99](https://github.com/PharmGKB/PharmCAT/commit/c247c99254427e22df4530ccdd656b49a1f4120b))
* **preprocessor:** flag non-PGx variants at PGx positions ([908f6d9](https://github.com/PharmGKB/PharmCAT/commit/908f6d9d83f867e7a496af0c25b73ce9f3cf29fb)), closes [#87](https://github.com/PharmGKB/PharmCAT/issues/87)


### Bug Fixes

* **data:** update to CPIC 1.14 ([e5ec432](https://github.com/PharmGKB/PharmCAT/commit/e5ec4324c8854e4c83817e995ac600e4b90d7e1c))
* **namedallelematcher:** add novel ALT warning ([0ed4462](https://github.com/PharmGKB/PharmCAT/commit/0ed44629527ad67dc0001791ac468ecceb453f8e))
* **namedallelematcher:** add support for PCATxALT and PCATxREF filters from preprocessor ([a81e35c](https://github.com/PharmGKB/PharmCAT/commit/a81e35c7e4ac5d5338a1f829c9ba1a9a67867c96))
* **preprocessor:** reading of input file list and sample list ([99d472c](https://github.com/PharmGKB/PharmCAT/commit/99d472c62c58fd92d1a065f8dee604bc2ceac9f7)), closes [#88](https://github.com/PharmGKB/PharmCAT/issues/88)
* **preprocessor:** remove key of chromosome positions only if the key is present in the dictionary ([ebf5ac2](https://github.com/PharmGKB/PharmCAT/commit/ebf5ac2d8d74d3472fe7c7530a4d045811fa6918))
* **preprocessor:** remove key of chromosome positions only if the key is present in the dictionary ([0177a38](https://github.com/PharmGKB/PharmCAT/commit/0177a384287a4ada2bbd97687b6453cd3124a0e6))
* **preprocessor:** remove key only if the key is present in the dictionary ([2f29fd1](https://github.com/PharmGKB/PharmCAT/commit/2f29fd1c31c4a54801a4681041242e08c2c47607))
* **preprocessor:** remove PASS if the position is flagged for PharmCAT; update FILTER flag descriptions ([f06c8c0](https://github.com/PharmGKB/PharmCAT/commit/f06c8c0900cbb84823da58f9f9ba513aea92268f))

## [1.4.0](https://github.com/PharmGKB/PharmCAT/compare/v1.3.1...v1.4.0) (2022-01-28)


### Features

* **data:** update data for clopidogrel guideline update ([b80a2ab](https://github.com/PharmGKB/PharmCAT/commit/b80a2ab745c9669a3939452f7a5f698f4b92e48d))
* **data:** update position data ([fad0485](https://github.com/PharmGKB/PharmCAT/commit/fad0485ee6218e57da143a2ebccc911b32ee979e))


### Bug Fixes

* **datamanager:** do not update docs if allele information not loaded ([e5ca7ac](https://github.com/PharmGKB/PharmCAT/commit/e5ca7accf7ee23677d5884a657712227bed34b64))
* **namedAlleleMatcher:** warn if ref in definition does not match ref in VCF ([e188e41](https://github.com/PharmGKB/PharmCAT/commit/e188e416bf8983558330402bcdc82e5c2b541419))
* **preprocessor:** detect VCF samples that contains ',' which violates the VCF convention ([59abc8a](https://github.com/PharmGKB/PharmCAT/commit/59abc8aa344e85703e7c7e7bfe1411479f029ef8))
* **preprocessor:** homozygous reference at a single-nucleotide locus will not infer the genotype status of the INDELs at the same genomic position ([261531b](https://github.com/PharmGKB/PharmCAT/commit/261531b24844df0eca7871644548be77a222c1f7))
* **preprocessor:** improve annotations of ID and info columns ([e0b8ba8](https://github.com/PharmGKB/PharmCAT/commit/e0b8ba87d8425e0a007e09b1c93e3defcc238040))
* **preprocessor:** remove redundant homozygous reference check ([bd8ba92](https://github.com/PharmGKB/PharmCAT/commit/bd8ba921140e9813a74febc86ef8eaa332e81aa0))

### [1.3.1](https://github.com/PharmGKB/PharmCAT/compare/v1.3.0...v1.3.1) (2022-01-14)


### Bug Fixes

* **data:** CPIC update ([097b275](https://github.com/PharmGKB/PharmCAT/commit/097b27597b2cb04aa83d1d0f58567c7eb0992845))

## [1.3.0](https://github.com/PharmGKB/PharmCAT/compare/v1.2.1...v1.3.0) (2021-12-09)


### Features

* **preprocessor:** improve support for phased data ("--phased") ([52f6ed0](https://github.com/PharmGKB/PharmCAT/commit/52f6ed0113445c9aef16e1972a71db85c92dc812)), closes [#75](https://github.com/PharmGKB/PharmCAT/issues/75) [#78](https://github.com/PharmGKB/PharmCAT/issues/78)


### Bug Fixes

* **data:** CPIC update ([e741c7f](https://github.com/PharmGKB/PharmCAT/commit/e741c7f991cbe0fda35839ceeef5340fc2112ceb))
* **preprocessor:** add missing multiallelic variants/positions as phased and bcftools determines phasing by GT delimiter accordingly ([f7762ad](https://github.com/PharmGKB/PharmCAT/commit/f7762adb6cb77b764dd2544183b7cb3b33e04f71)), closes [#78](https://github.com/PharmGKB/PharmCAT/issues/78)
* **preprocessor:** fix for PGx positions with missing ALT ([e2592b2](https://github.com/PharmGKB/PharmCAT/commit/e2592b2e3f12a44566140753aac51a8d8bcc9240)), closes [#77](https://github.com/PharmGKB/PharmCAT/issues/77)
* **preprocessor:** fix output dir of a temp ([f726303](https://github.com/PharmGKB/PharmCAT/commit/f72630345db80920a8869127564e6cf87b628f98))
* **preprocessor:** interrupt and print a warning message if a gVCF input is detected ([441a384](https://github.com/PharmGKB/PharmCAT/commit/441a38474b60e1918e4a7ce0ca9d7198b5408b02))

### [1.2.1](https://github.com/PharmGKB/PharmCAT/compare/v1.2.0...v1.2.1) (2021-11-18)


### Bug Fixes

* **data:** fix bad genotype in a CYP2C19 test case ([0a77c77](https://github.com/PharmGKB/PharmCAT/commit/0a77c7715ee86a81bbccebd202ae488ec83e31a2))
* **DataManager:** support wobble code in reference allele ([ada5a3a](https://github.com/PharmGKB/PharmCAT/commit/ada5a3a67e091ba53a309966144722c2f520c995))
* **DataManager:** treat HGVS dup as a form of repeat [2] ([81135db](https://github.com/PharmGKB/PharmCAT/commit/81135db70719a0f1e71e62e86e9687f0c550c7ca))
* **preprocessor:** fix vcf header parsing error ([b5c9814](https://github.com/PharmGKB/PharmCAT/commit/b5c981425d12856c15bf6ea39145e276f5a73ca5))

## [1.2.0](https://github.com/PharmGKB/PharmCAT/compare/v1.1.0...v1.2.0) (2021-10-27)


### Features

* set missing positions to ref ([d80e8a7](https://github.com/PharmGKB/PharmCAT/commit/d80e8a7c0be39c5ccc83b5d9b85f4e004ea52d39))
* **data:** update to CPIC v1.10 ([0547340](https://github.com/PharmGKB/PharmCAT/commit/05473404505229ce76493552a56615ff5d3a017a))


### Bug Fixes

* remove UGT1A1 special handling and add more SLCO1B1 examples ([74c9ae8](https://github.com/PharmGKB/PharmCAT/commit/74c9ae8004ff20c6153ec3d170c373ccc43864ff))
* update vcf header and properly sort vcf after normalization ([e512bc3](https://github.com/PharmGKB/PharmCAT/commit/e512bc321c86ae14596f4475bf5423f757bb5e24))
* validate bgzip ([78ccd71](https://github.com/PharmGKB/PharmCAT/commit/78ccd7187f12cd7c875a358da84e1efe8b525bf1))
* **data:** update SLC01B1 data ([47b572a](https://github.com/PharmGKB/PharmCAT/commit/47b572a3ebdb4e2bb80e0a23d145a6a1d064d33c))
* **docker:** update docker to use bgzipped reference FASTA ([15e5977](https://github.com/PharmGKB/PharmCAT/commit/15e597796cfc115f7681de2c8124890c79855ac3))
* **docker:** update docker to use bgzipped reference FASTA ([9b594f6](https://github.com/PharmGKB/PharmCAT/commit/9b594f6b6b3e1ecd3d13df88c70986851a51e1fc))
* **NamedAlleleMatcher:** improve warning message when GT in VCF doesn't have 2 alleles ([890ae92](https://github.com/PharmGKB/PharmCAT/commit/890ae9206d25a527844d7a65f55e19ca1c9d4043))
* **preprocessor:** if output_folder is not specified, use parent directory of input ([fef2722](https://github.com/PharmGKB/PharmCAT/commit/fef2722260ee5a15151a2345b1a44e1b4c6f039f))
* **preprocessor:** improve how reference FASTA is obtained ([9860ce2](https://github.com/PharmGKB/PharmCAT/commit/9860ce28df9844058a20bcad9bf21d2a3c9d4fbb))
* **preprocessor:** update usage docs ([0faf08e](https://github.com/PharmGKB/PharmCAT/commit/0faf08ea8126b414e8995a2c3da2e90f4cc2628c))

## [1.1.0](https://github.com/PharmGKB/PharmCAT/compare/v1.0.0...v1.1.0) (2021-10-14)


### Features

* default '--ref_pgx_vcf' to 'pharmcat_positions.vcf.bgz' in the current working directory ([eb34576](https://github.com/PharmGKB/PharmCAT/commit/eb3457642abd6a74ba1e2596624398cfb23abf2e))
* make output directory optional ([55c244b](https://github.com/PharmGKB/PharmCAT/commit/55c244b62564fffb01ca9d4da53862bd744f364b)), closes [#68](https://github.com/PharmGKB/PharmCAT/issues/68)
* normalize "chrMT" to "chrM" ([509c010](https://github.com/PharmGKB/PharmCAT/commit/509c010bf8381d46cda1e614f1cff4feda3e67e6))
* output dir of preprocessor now defaults to the dir of input VCF ([c677e61](https://github.com/PharmGKB/PharmCAT/commit/c677e6175e7e4aaa1847cd0f7334d24dc849869a))


### Bug Fixes

* display version with --version flag ([f7f6bbe](https://github.com/PharmGKB/PharmCAT/commit/f7f6bbe9fe3f7808dba8e6a9863fe6c56d395ac4))
* include URL to docs when multisample VCF is found ([6105bcb](https://github.com/PharmGKB/PharmCAT/commit/6105bcb836124aabd9d52cd5a50cb57816da4551))
* normalize "chrom" field to chrM for mitochondria ([da2a456](https://github.com/PharmGKB/PharmCAT/commit/da2a456cd2b7707f55b3c91af99698ccca3c9676))
* remove pre-release note in PharmCAT runtime ([81e9233](https://github.com/PharmGKB/PharmCAT/commit/81e9233b821baebbcd619f04c2b21f7e3b6b1fe7))
* sort pharmcat_positions.vcf, add support for .bgz and .tbi of pharmcat_positions.vcf ([2434f56](https://github.com/PharmGKB/PharmCAT/commit/2434f5634cb1ed1d45484e720e103a254024f248))
* update processor details to v1.0.0 ([ac3d001](https://github.com/PharmGKB/PharmCAT/commit/ac3d001fe35b6fa01064a0877e92ed4bc9b1f385))

## [1.0.0](https://github.com/PharmGKB/PharmCAT/compare/v0.8.0...v1.0.0) (2021-09-27)


### ⚠ BREAKING CHANGES

* The allele definition format has been updated and is not backwards compatible.
* updating to Java 14

### Features

* add "pj" flag to PharmCAT class for writing Phenotyper output to a JSON file ([5986b66](https://github.com/PharmGKB/PharmCAT/commit/5986b66c9af7b1543df34022a4962bb933e39516))
* add a summary page to PharmCAT website ([565f031](https://github.com/PharmGKB/PharmCAT/commit/565f0313f39d693b85a057aa5ce67210f7efbf5f))
* add AutogeneratedVcfTester ([e504c24](https://github.com/PharmGKB/PharmCAT/commit/e504c24871a175f56dd3f85bf2385cad16c04030))
* add CLI option to get all results from NamedAlleleMatcher ([40c4377](https://github.com/PharmGKB/PharmCAT/commit/40c437797151e7b539cec020778e550ae21e0232))
* add CPIC version to Reporter output ([83027f4](https://github.com/PharmGKB/PharmCAT/commit/83027f4c58c189a105be5c33e90f2c7d7401a72d))
* add docker support ([b53f46d](https://github.com/PharmGKB/PharmCAT/commit/b53f46de5849feb6016db8ef00274180c3effae9))
* add drugs file to DataManager and update drugs file ([d604652](https://github.com/PharmGKB/PharmCAT/commit/d604652f819a6d192a9fd844563c1e920926ae6a))
* add exact-match-only to AutogeneratedVcfTester ([fd14cd6](https://github.com/PharmGKB/PharmCAT/commit/fd14cd6023fa0269a94e3c0b074ca6f1beb50c62))
* add messages to drug reports for certain gene calls ([3da72a0](https://github.com/PharmGKB/PharmCAT/commit/3da72a0fe5c5d5640c877482f71425d18aefe192))
* add option to PharmCAT runner to retain all scoring matches ([f8b8ece](https://github.com/PharmGKB/PharmCAT/commit/f8b8ece124fa5e687ea932cd672a4ed922ab7765))
* add support for MT-RNR1 ([6d93598](https://github.com/PharmGKB/PharmCAT/commit/6d93598d26df1c8c6b12e74a3fdcda4939a4cd32))
* add warning messages for ambiguous diplotype calls ([2c2cc37](https://github.com/PharmGKB/PharmCAT/commit/2c2cc37456b7ed9149970075b92418356c4b51da))
* data update from CPIC and related changes ([3d23a85](https://github.com/PharmGKB/PharmCAT/commit/3d23a857b5f9d3cc24b8d21f6dff72c5bac36cee))
* extract only exactly matched PGx variants (used to based on position) ([e101d27](https://github.com/PharmGKB/PharmCAT/commit/e101d278a11344f5e714562819121c3c5de40af2))
* handle unassigned function alleles ([cebfda9](https://github.com/PharmGKB/PharmCAT/commit/cebfda9d5b605aca218cc2953a4bcb3caacf4ce2))
* make an explicit list of "reportable" drugs ([99bf3b4](https://github.com/PharmGKB/PharmCAT/commit/99bf3b45fae51606034071707881383cac457821))
* move CYP2D6 to list of preferred outside calls in summary report ([8d07f60](https://github.com/PharmGKB/PharmCAT/commit/8d07f60bec3a8fbe5eec812cd38b7f55030fb3f5))
* new unphased data note in final report ([826cbb5](https://github.com/PharmGKB/PharmCAT/commit/826cbb5752860e5e11239f5857a06a4c8327f484))
* normalize alleles for VCF ([bfd37d8](https://github.com/PharmGKB/PharmCAT/commit/bfd37d8614f8a2bd5705c99a2cc1a0f7be00589d))
* option to keep intermediate files ([882990b](https://github.com/PharmGKB/PharmCAT/commit/882990b3ce494cc39232586790ddf80dd1a6d331))
* option to provide a list of vcf files ([dbef0c4](https://github.com/PharmGKB/PharmCAT/commit/dbef0c4d9e95e86cac0d5b6d76d27121624948ab))
* show matching diplotypes in final report recommendations sections ([de18261](https://github.com/PharmGKB/PharmCAT/commit/de18261d117c8931b1a897be8c22509f1c50772f))
* support new message annotation matching for ambiguous het calls ([831fb79](https://github.com/PharmGKB/PharmCAT/commit/831fb798822d042867ac5188f8e72b64e2a01321))
* take a file of samples to preprocess ([b0f50f2](https://github.com/PharmGKB/PharmCAT/commit/b0f50f26469d295e839337573a68488ed3d76015))
* update report disclaimer template to match website ([75d1af3](https://github.com/PharmGKB/PharmCAT/commit/75d1af321e4133870e48669e7738e604879ab3ae))
* validate bcftools and tabix ([252f7c1](https://github.com/PharmGKB/PharmCAT/commit/252f7c11237d15e1c7d185be8a6cea12065de4dd))


### Bug Fixes

* add fuzzy match support to AutogeneratedVcfTester ([296fef6](https://github.com/PharmGKB/PharmCAT/commit/296fef6a2a06ea94a940be6245090ecd13754b68))
* add overlap check to AutogeneratedVcfTester ([bad7b98](https://github.com/PharmGKB/PharmCAT/commit/bad7b986c348d9567bcb4d3ef2078c6b7c6b2f3b))
* ambiguous code expansion ([f097ce7](https://github.com/PharmGKB/PharmCAT/commit/f097ce76877f62956e2a350c39f3e0d9ffefd75e))
* ambiguous REF in output VCF files ([d8a9c17](https://github.com/PharmGKB/PharmCAT/commit/d8a9c1739730ac26b278cb347ead19a6ae025090))
* apply same criteria for ambiguity messages to genes and drugs ([fa94e2a](https://github.com/PharmGKB/PharmCAT/commit/fa94e2a37fd8abac4ab4893d1f8d07bcf53355ea))
* change naming convention for IFNL3 to IFNL3/4 for final report ([ffd1aa8](https://github.com/PharmGKB/PharmCAT/commit/ffd1aa8507cd02eeafc527f79e49e16bd0430fe2))
* clean up message annotation matching ([1c0b4f4](https://github.com/PharmGKB/PharmCAT/commit/1c0b4f4c72e323f16204873cf3c4eb3bb66c9316))
* cleanup, appease linter ([b260c24](https://github.com/PharmGKB/PharmCAT/commit/b260c24c4a9d33b0b174db95ab306e29a70499e2))
* custom definition transform for CYP2C19 *1 and *38 ([1489005](https://github.com/PharmGKB/PharmCAT/commit/148900585d34038937c427e77fabcb77ae8ded9e))
* default "No Result" value for uncalled genes when doing recommendation lookup ([60961a8](https://github.com/PharmGKB/PharmCAT/commit/60961a80e1e7f27d48bbd628cab1cf5a0315f852))
* default allHits/assumeReference in DefinitionExemption to be null ([fcc23f2](https://github.com/PharmGKB/PharmCAT/commit/fcc23f2f0c3b731dcb347650d651c04c2f0f884f))
* download from url, test vcf ([5840864](https://github.com/PharmGKB/PharmCAT/commit/58408646e9d9f247cf1aac1c289da3bc81e7a6d7))
* error in reading the sample file ([9dccb85](https://github.com/PharmGKB/PharmCAT/commit/9dccb852fbf4ba9f71832a7ace724b3f74aeffe9))
* exclude CYP2D6 from ExtractPositions ([71ef0cb](https://github.com/PharmGKB/PharmCAT/commit/71ef0cbb7cd52df52168d7be278d4b165e4139fe))
* file path split ([62427e9](https://github.com/PharmGKB/PharmCAT/commit/62427e9420a0ffd903cc8481a4b9a47a99aba08b))
* files containing wobble correctly identified ([a561daf](https://github.com/PharmGKB/PharmCAT/commit/a561dafc5f2d9d0ffda534863d534f89d811c9ab))
* finalize vcf query caching ([b031691](https://github.com/PharmGKB/PharmCAT/commit/b0316912451a38b745ac917344ad49da38250a77))
* fix ambiguity criteria for message annotation ([2b0074e](https://github.com/PharmGKB/PharmCAT/commit/2b0074eb2204bc6f2033523bb86b387934be890c))
* fix bad import ([4642716](https://github.com/PharmGKB/PharmCAT/commit/46427163b40ee1ddea7628abfa41b1882344042b))
* fix COPY error in Dockerfile ([37ddbfa](https://github.com/PharmGKB/PharmCAT/commit/37ddbfa00e4b32f567bd6fe727955d2295ab52a6))
* fix executable flag for test script ([5274440](https://github.com/PharmGKB/PharmCAT/commit/5274440a57634d8353d191b6c80be17e0397e1e0))
* fix ExtractPositionsTest to also ignore genes ([60d7e90](https://github.com/PharmGKB/PharmCAT/commit/60d7e906bc09f2ce3d71d92095e5d343ccadd199))
* fix format of examples page on site ([f8f4216](https://github.com/PharmGKB/PharmCAT/commit/f8f4216994c2872da36f3764fc2ebe5b5984a887))
* fix gradle "dataUpdate" task file paths to currently used paths ([90e7372](https://github.com/PharmGKB/PharmCAT/commit/90e7372d74194676f2be2b2c3424553d99e5e7ff))
* fix missing variant column in genotype table of report ([512742a](https://github.com/PharmGKB/PharmCAT/commit/512742ae4a5220c0350572ee9972a926c2e11ee5))
* fix NPE when variant has no dbSNP ID ([047313d](https://github.com/PharmGKB/PharmCAT/commit/047313d6085af2bb70371cb0531efc87657b8612))
* fix pre-release URL ([a63e703](https://github.com/PharmGKB/PharmCAT/commit/a63e703d58a7c8f568f4f23e6438a26c1c320583))
* fix report layout ([59b8e1c](https://github.com/PharmGKB/PharmCAT/commit/59b8e1c55e2ea5af55443708f4363bac96277ebd))
* fix sequence for deletion in NUDT15 definition ([9a369d7](https://github.com/PharmGKB/PharmCAT/commit/9a369d78a06d59c40e4cc1a55e5728d4c3b2acda))
* fix summary report of "outside call" genes and update summary page ([a8e341c](https://github.com/PharmGKB/PharmCAT/commit/a8e341ca106b5c8a3bb213c538c751c217728034))
* fix the image URL for warfarin diagram ([643317c](https://github.com/PharmGKB/PharmCAT/commit/643317c4b2b001d412223034e175b5c25c80f7af))
* fix typo in message annotation data file ([3b93329](https://github.com/PharmGKB/PharmCAT/commit/3b9332939638abb3e0be25e1b26cec164e3549bf))
* fix warfarin display of recommendation text ([57eaf10](https://github.com/PharmGKB/PharmCAT/commit/57eaf10bc02bef4ef5eae4d03338ab953d4f7659))
* fix wording for outside call alert ([61419c9](https://github.com/PharmGKB/PharmCAT/commit/61419c990127a15f14bfdafc2cd31094fd4c4076))
* fold pharmcat_positions.vcf into DataManager, remove pharmcat_intervals.txt ([bf0caa5](https://github.com/PharmGKB/PharmCAT/commit/bf0caa5aa5ad23bcb7ab63077a6eaedd232954f4))
* generate uncompressed vcf ([91bd1d2](https://github.com/PharmGKB/PharmCAT/commit/91bd1d2af1e864ce7c93ac6149fa47395657bb9b))
* handle ambiguous IUPAC codes properly ([fef7251](https://github.com/PharmGKB/PharmCAT/commit/fef72519ffdd1b1e2d32ff5b5a2d3e3ee08b3b0d))
* handle long list of biobank samples ([078f688](https://github.com/PharmGKB/PharmCAT/commit/078f6880cccc5510c89950b3aeac75a53b48d98c))
* handle unexpected alleles gracefully ([2dfbf0e](https://github.com/PharmGKB/PharmCAT/commit/2dfbf0ec817bf07538d3f11489460bceeb1e54ca))
* ignore MT-RNR1 on the reporter side ([ed4ddb1](https://github.com/PharmGKB/PharmCAT/commit/ed4ddb14f9c3784409ca9ed1b7ae5f539ecb606b))
* ignore positions with empty reference allele ([b5198b0](https://github.com/PharmGKB/PharmCAT/commit/b5198b0817c87645957c865cd2b97a226c799674))
* import package ([69ea88e](https://github.com/PharmGKB/PharmCAT/commit/69ea88e91edaae194871400d35d8a1ea80056c81))
* improve error message on reference allele mismatch, update positions_reference.tsv ([e4f414d](https://github.com/PharmGKB/PharmCAT/commit/e4f414d032d3595b12540e54a3b7c5fd72b4e955))
* indel ALT/REF and positions in output VCF ([2f311e5](https://github.com/PharmGKB/PharmCAT/commit/2f311e594b31f779ab240988214783babb29a8bd))
* make gene section name consistent in final report ([f8453ce](https://github.com/PharmGKB/PharmCAT/commit/f8453ce8cdf83787964cfac6e0f26421adfb616a))
* make het function phrases display consistently ([d800c69](https://github.com/PharmGKB/PharmCAT/commit/d800c69567b2fa15fc3ee7921b4eaf3ddc092606))
* make Reporter follow the modular pattern ([4b4f567](https://github.com/PharmGKB/PharmCAT/commit/4b4f5679d1ef2c49618625f2964afbaf335342ae))
* merge vcf position extraction and intervals into DataManager ([ee7e70e](https://github.com/PharmGKB/PharmCAT/commit/ee7e70eb8730631d5213308a9042901a2012cf6c))
* modify tabix error message ([9a0fd2d](https://github.com/PharmGKB/PharmCAT/commit/9a0fd2d970739a3fb4d1e7be3e7fe8d4ebe9153a))
* preprocessing test file format errors ([b4234b1](https://github.com/PharmGKB/PharmCAT/commit/b4234b1dbc234765f855c9764b6e6581538f5b78))
* prevent *1 from being shown as a relevent allele for all positions in CYP2C19 ([57b8b27](https://github.com/PharmGKB/PharmCAT/commit/57b8b27833a2593462e33b8a262a7a97325e2e94))
* print out the full error traceback message ([a39b6c5](https://github.com/PharmGKB/PharmCAT/commit/a39b6c55820cb966e3f0e19667763d11465c0d27))
* pull Sheets by URL instead of using Sheets API ([5c1b482](https://github.com/PharmGKB/PharmCAT/commit/5c1b4823923e7ac9315de783199c85a62bc191e5))
* remove "show-all-matches" option from PharmCAT class ([062edec](https://github.com/PharmGKB/PharmCAT/commit/062edeccdd4561bc5c450440fc2f8f543bb23915))
* remove cyvcf2 package that is not compatible with FIPS ([37ec7e9](https://github.com/PharmGKB/PharmCAT/commit/37ec7e972f2ce66a888acb0f100f12321c7e4ecc))
* remove duplicate recommendations from Reporter output ([0d319c5](https://github.com/PharmGKB/PharmCAT/commit/0d319c5231ff934bf2d05f5d8aaffd1a8bd111ce))
* remove references to unneeded pharmcat.properties file ([66d726e](https://github.com/PharmGKB/PharmCAT/commit/66d726efb45948f35072d85b44c87a5f1e1d6cd8))
* remove the import of cyvcf2 ([62fc70c](https://github.com/PharmGKB/PharmCAT/commit/62fc70cf9a18935bbaa3c42f79c73c3773c39881))
* remove unnecessary footnote from genotypes table in final report ([191e469](https://github.com/PharmGKB/PharmCAT/commit/191e469a82f20c934c3c5a19f4bfa766b14fb93c))
* restore default "assumeReference" behavior to UGT1A1 ([8e871bb](https://github.com/PharmGKB/PharmCAT/commit/8e871bbd8b625b6a26a9d94b064b397fe5d92d37))
* show DiplotypeMatch names in order of haplotype names ([83e31e5](https://github.com/PharmGKB/PharmCAT/commit/83e31e5b89cdcdcca4413c6f4cb57451449d797c))
* show unphased footnote superscript on applicable calls ([d28e181](https://github.com/PharmGKB/PharmCAT/commit/d28e1819b9e157bfada0b4569677c451e634ffd4))
* support "." in VCF allele field ([fbd20f0](https://github.com/PharmGKB/PharmCAT/commit/fbd20f0d3381777b9347857abd914a422628e254))
* support renaming downloaded fasta files ([1c5ff92](https://github.com/PharmGKB/PharmCAT/commit/1c5ff92de8713600cb732225e75e656cee7a4c7d))
* support single allele call in VCF ([988771b](https://github.com/PharmGKB/PharmCAT/commit/988771b0e145fc68328dd38ddfd672d6f47f6a25))
* switch named allele collection from List to SortedSet ([c734716](https://github.com/PharmGKB/PharmCAT/commit/c7347162e0918a378229b2ae74971d3a3ff2881b))
* take sorted by-chromosome VCFs as input ([18a7fc9](https://github.com/PharmGKB/PharmCAT/commit/18a7fc920a7c36cd472624d57228f37747fef119))
* test_gen scripts no longer need positions.vcf ([d593d6e](https://github.com/PharmGKB/PharmCAT/commit/d593d6e6e675cfd1e790ef9119d1f20c24210b9a))
* update and fixes for allele definition data ([ce11f6d](https://github.com/PharmGKB/PharmCAT/commit/ce11f6ddad72a1424cfc2ad42cc48053d436ac60))
* update copy for missing alleles in gene match section ([4540cad](https://github.com/PharmGKB/PharmCAT/commit/4540cad002d5c382004dc4e435516bd2841a3b14))
* update CYP2C19 to account for new reference haplotype ([e9856e1](https://github.com/PharmGKB/PharmCAT/commit/e9856e140485a95f7459fdaac53ae736ee6f2330))
* update definitions ([47c549f](https://github.com/PharmGKB/PharmCAT/commit/47c549f6f8878ef0152be0b5e09d91f9aca771e0))
* update definitions/exemptions/messages ([9a45d05](https://github.com/PharmGKB/PharmCAT/commit/9a45d0513c4848c090566181cd7008ca9c9f8cfd))
* update shebang to use python3 ([f103fd6](https://github.com/PharmGKB/PharmCAT/commit/f103fd68ffc60343e3a9d182a648603d41c383cf))
* VCF filename matching with wobbles across runs ([0b77cde](https://github.com/PharmGKB/PharmCAT/commit/0b77cdeba64525c4c3e0eb9e300b285653a1aa05))


### Performance Improvements

* optimize log messages ([a79ec99](https://github.com/PharmGKB/PharmCAT/commit/a79ec99b4e6877133693445ab8522639bcd58744))


### Miscellaneous Chores

* updating to Java 14 ([c91edeb](https://github.com/PharmGKB/PharmCAT/commit/c91edeb103a5d056caa61273f773cfb140ab9ea9))

## [0.8.0](https://github.com/PharmGKB/PharmCAT/compare/v0.7.0...v0.8.0) (2021-03-31)


### Features

* add "label" field to Phenotyper output for diplotypes ([b92b914](https://github.com/PharmGKB/PharmCAT/commit/b92b9144ba8b90c62c143c912a69283a44c9529e))
* add new GenotypeInterpretation class ([c90e358](https://github.com/PharmGKB/PharmCAT/commit/c90e3580f8562130e547087fc01c99801aed20e8))
* add the reference flag to namedAllele objects ([f02d872](https://github.com/PharmGKB/PharmCAT/commit/f02d87287460fe61db05746e84fe080e3d37f1b5))
* add timer ([eb5c997](https://github.com/PharmGKB/PharmCAT/commit/eb5c99790a9e759c8f8abee788eedc2add7a930b))
* add warfarin and peginterferon back in to drug list ([18ed7f9](https://github.com/PharmGKB/PharmCAT/commit/18ed7f94334265e02fdb1bdd5c88e1bf349e1bba))
* change gene-phenotype to use new layout ([151de09](https://github.com/PharmGKB/PharmCAT/commit/151de0949add7f9a1846673e42e44cd503348d49))
* download to temp file ([5c46985](https://github.com/PharmGKB/PharmCAT/commit/5c46985438deddd3fdd49cd58e69c4f5ee8975f7))
* enable Phenotyper to take either VCF input or NamedAlleleMatcher output JSON ([1631d52](https://github.com/PharmGKB/PharmCAT/commit/1631d528f0f494eadf450e537359576e7e6a4acf)), closes [#39](https://github.com/PharmGKB/PharmCAT/issues/39)
* handle local files ([9d0c46b](https://github.com/PharmGKB/PharmCAT/commit/9d0c46b33c89c355c0ec7fcaecd74bba6fb6526c))
* improve and expand outside call features ([e2eedfc](https://github.com/PharmGKB/PharmCAT/commit/e2eedfcbfa81a2bea7f6084a4fb14f76400c87b6))
* include "comments" on recommendations in final report ([71ce411](https://github.com/PharmGKB/PharmCAT/commit/71ce411a6bcb57da9d6d8c42890892c2449f8fb2))
* pull definitions from S3 and support ignored positions ([4702f04](https://github.com/PharmGKB/PharmCAT/commit/4702f04e4c86e0bb65a7633ca26b2b5b7b36cfa6))
* recursive output directory creation ([b547294](https://github.com/PharmGKB/PharmCAT/commit/b54729441ababa4d0b10b3cad7d0a65d5850db50))
* reporter update intermediate check-in ([33c560b](https://github.com/PharmGKB/PharmCAT/commit/33c560bdad8797fea5100fbd136962a5784ddbe2))
* update gene definition files ([141bed6](https://github.com/PharmGKB/PharmCAT/commit/141bed6125654b4868f1658f177a110d3c26ade8))
* update to latest CPIC drugs and phenotypes ([780a1f6](https://github.com/PharmGKB/PharmCAT/commit/780a1f68b4dde093a3bbbf829cd237116813ca8e))


### Bug Fixes

* adjust TPMT variant count ([dceb6ba](https://github.com/PharmGKB/PharmCAT/commit/dceb6ba28f6d8da344054e2dc9b0723dd50860d4))
* appease javadoc ([6153145](https://github.com/PharmGKB/PharmCAT/commit/615314546aa2f4e4d3be24150ed52cb96d78735d))
* avoid warning for transient Pattern in NamedAllele ([4e27107](https://github.com/PharmGKB/PharmCAT/commit/4e2710716a5fd06f334507a8526e4b3c272ce8a3))
* clear rsid map before rebuilding it ([d9973b8](https://github.com/PharmGKB/PharmCAT/commit/d9973b8b45eaefed88e85d3b06067a150eb4d865))
* combine duplicate position columns from CYP2D6 translation ([e70648d](https://github.com/PharmGKB/PharmCAT/commit/e70648d9285c4006a620711fe9898911ef93cd13))
* display number of files produced ([842414d](https://github.com/PharmGKB/PharmCAT/commit/842414db9ca5a8dcc2f37e66a6f0a1c79d81b378))
* do not include data for anything involving blacklisted genes ([408f503](https://github.com/PharmGKB/PharmCAT/commit/408f503be7bf1eebbf1b31c245cbd5382490df0a))
* DPYD chromosomal HGVS name ([d82cb44](https://github.com/PharmGKB/PharmCAT/commit/d82cb44df8efa2722479b433f85934f1c23d5510))
* explicitly set date format ([876182e](https://github.com/PharmGKB/PharmCAT/commit/876182e6c1914f686285f6764855343ffb24c4a7))
* fix bug when diplotypes are specified in reverse order (fixes DPYD bug) ([3849c47](https://github.com/PharmGKB/PharmCAT/commit/3849c478c79e2a34793c591f9e52dbcba814c40c))
* fix bug when overriding callable genes with outside calls ([85ad9f5](https://github.com/PharmGKB/PharmCAT/commit/85ad9f5de54dbdbbd3a9198ab5fc38e5ef4deae8))
* fix bugs and performance in ExtractPositions ([95f3725](https://github.com/PharmGKB/PharmCAT/commit/95f3725dba91afe34ce142338777090af4315fc4)), closes [#34](https://github.com/PharmGKB/PharmCAT/issues/34)
* fix CFTR tests to remove F508del and change "Reference" ([c50cfd4](https://github.com/PharmGKB/PharmCAT/commit/c50cfd4d1d5a79227b2c10b974a35f99972e3fa5))
* fix CYP2C9 *2/*3 unit test ([d82fd7e](https://github.com/PharmGKB/PharmCAT/commit/d82fd7e1b1fcba56eb12b4d4a26ce6f5102428fd))
* fix delete obsolete files ([11da562](https://github.com/PharmGKB/PharmCAT/commit/11da5620e244df6af718d8651b0bda016fcf01fd))
* fix DPYD integration test ([524c3c3](https://github.com/PharmGKB/PharmCAT/commit/524c3c3865818b81fc9f2f4ae714827ca8259dcb))
* fix gene definition files to make reference allele be first ([0951f57](https://github.com/PharmGKB/PharmCAT/commit/0951f579957e9f2cf24f8f64f4316c42bc235279))
* fix NamedAlleleMather tests to account for changed allele names ([14f0f7f](https://github.com/PharmGKB/PharmCAT/commit/14f0f7f497a569cc1502d56718704232025a07bf))
* fix recommendation matching for multi-gene guidelines ([bd1ee0c](https://github.com/PharmGKB/PharmCAT/commit/bd1ee0c0cfdf8e77a4bc725dd1547d819f543eeb))
* fix remove ignored positions not removing associated alleles ([9ca4a2c](https://github.com/PharmGKB/PharmCAT/commit/9ca4a2cc34c9f1d5b804b1986b816a6a2dcfaee0)), closes [#36](https://github.com/PharmGKB/PharmCAT/issues/36)
* fix unsafe operation compiler warning ([a081bb7](https://github.com/PharmGKB/PharmCAT/commit/a081bb74985248138da133f50e24a6c40331f2e2))
* fix variant ordering bug ([ee7e1af](https://github.com/PharmGKB/PharmCAT/commit/ee7e1af0b3f697e92c75de06499281454ea471be))
* fix VCF syntax for unspecified and deletion genotypes ([5b6ae8e](https://github.com/PharmGKB/PharmCAT/commit/5b6ae8ec015327a8064e8d33e8bb182016167841))
* improve handling of "Unknown" calls ([a46866f](https://github.com/PharmGKB/PharmCAT/commit/a46866f894fa2a02a779440b4e88f9b0c44b88d2))
* names for position chr6:18132163 in TPMT ([39a99f7](https://github.com/PharmGKB/PharmCAT/commit/39a99f79b5773c6e4fedb518d723beb90dd6fe6f))
* remove "highlighted" drugs in final report ([4a1f9ef](https://github.com/PharmGKB/PharmCAT/commit/4a1f9efb8a47bf48c049c61bce480a0d507989ef))
* remove "no calls" from sample files ([1fbb002](https://github.com/PharmGKB/PharmCAT/commit/1fbb002379f6631b0d153ef29d81caaa9ac64dd9))
* remove genes not used in recommendation match ([7c027d9](https://github.com/PharmGKB/PharmCAT/commit/7c027d91cf864d1e65407776d43e0ef642d3c55e))
* remove ignored positions from named alleles as well ([351979f](https://github.com/PharmGKB/PharmCAT/commit/351979f2c56ed23e2db73ddb2217c97041e8b8bb))
* remove styling for Rx change in report ([4e89125](https://github.com/PharmGKB/PharmCAT/commit/4e891252dc6e7b9c20b28b9cb5a65d23ae74a2f3))
* remove tests for *60 UGT1A1 allele ([80fbd5f](https://github.com/PharmGKB/PharmCAT/commit/80fbd5f18795f14a1f5a549432aa2b5414297c0b))
* remove unnecessary variant allele options ([f931b89](https://github.com/PharmGKB/PharmCAT/commit/f931b89a02e31a55369e7acae2677c8855cffcc0))
* remove unused "g" option for PharmCAT CLI ([f161a67](https://github.com/PharmGKB/PharmCAT/commit/f161a672e2db442bc32dec7bee2e25da159fb988))
* removed redundant convert_to_*.py scripts; condensed the scripts to one main script and one function library ([8ef9f21](https://github.com/PharmGKB/PharmCAT/commit/8ef9f2189a86a8075ec4cb3362e5b292f442b315))
* show genotype of "highlighted" variants in guideline section of report ([f2bf726](https://github.com/PharmGKB/PharmCAT/commit/f2bf726ef06f21b39f67ab18aed52c6195e88113)), closes [#31](https://github.com/PharmGKB/PharmCAT/issues/31)
* skip import of allele definitions that have all alleles ignored ([65dacd9](https://github.com/PharmGKB/PharmCAT/commit/65dacd9f660673e55a8ad92d342179a172f23df5))
* support ignored position exemptions ([5ecb6e1](https://github.com/PharmGKB/PharmCAT/commit/5ecb6e1b5991931146f12ee93a5d076f424af068))
* switch to PEP8 code style, clean up output, add basic input error checking ([dac626c](https://github.com/PharmGKB/PharmCAT/commit/dac626ca7272a064e075e805949782bb6625cd05))
* take chr## or ## (-> chr##) for CHROM ([9bce353](https://github.com/PharmGKB/PharmCAT/commit/9bce35313ef802f3dfc1a8174699ced9f4803478))
* update input argument ([2f57352](https://github.com/PharmGKB/PharmCAT/commit/2f573526a92befba7e23ba22ae6ca71b5978dace))
* update test cases to adjust for changes to allele definitions and exemptions ([cc0d997](https://github.com/PharmGKB/PharmCAT/commit/cc0d99733cb0fa62b1aaf0aa90796b4ab6d99dcd))


### Performance Improvements

* speed up VCF preprocessing ([d5a847d](https://github.com/PharmGKB/PharmCAT/commit/d5a847d147529482301197974407ddd5fc740414))
