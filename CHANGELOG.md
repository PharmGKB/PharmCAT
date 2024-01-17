## [2.9.0](https://github.com/PharmGKB/PharmCAT/compare/v2.8.3...v2.9.0) (2024-1-17)


### Features

* **data:** update to CPIC version v1.35 ([5cb6717](https://github.com/PharmGKB/PharmCAT/commit/5cb671756393b82732c2ffeb5e06201d453cade3))
* **data:** update to PharmGKB version 2023-12-19 ([66c49f9](https://github.com/PharmGKB/PharmCAT/commit/66c49f9d584b9c3e72cb6488fc17ffc79e9442cd))
* **reporter:** add new "non-match" type of message ([f5df05a](https://github.com/PharmGKB/PharmCAT/commit/f5df05a0cf0c7edd8e91fc59644cb81ea0535504))


### Bug Fixes

* **namedallelematcher:** DPYD matching should ignore find-combinations mode ([153853a](https://github.com/PharmGKB/PharmCAT/commit/153853a36a79694002de3d92c5fd56ad17a113d2))
* **namedallelematcher:** improve AD number warning ([9d9a0b9](https://github.com/PharmGKB/PharmCAT/commit/9d9a0b98dbd5ea75ea5ab7476986f715125eb0a3)), closes [#168](https://github.com/PharmGKB/PharmCAT/issues/168)
* **phenotyper:** fix bug with stripping gene symbol from outside call diplotypes ([cb97898](https://github.com/PharmGKB/PharmCAT/commit/cb978987fc7e54386c1bd09e6caeba817e485de2)), closes [#161](https://github.com/PharmGKB/PharmCAT/issues/161)


### Performance Improvements

* improve runtime of diplotype comparison ([406930c](https://github.com/PharmGKB/PharmCAT/commit/406930cedffa551ed3e21331766c8f630a9b7881))

## [2.8.3](https://github.com/PharmGKB/PharmCAT/compare/v2.8.2...v2.8.3) (2023-10-24)


### Bug Fixes

* improve labeling of reference calls ([4713b9b](https://github.com/PharmGKB/PharmCAT/commit/4713b9bae11e2de782e682bc2a8c935a12b679a2)), closes [#158](https://github.com/PharmGKB/PharmCAT/issues/158)
* **data:** fix date parsing for Java 21 ([aba3dcd](https://github.com/PharmGKB/PharmCAT/commit/aba3dcd8cdd0377175038a51c1435ac523141fe3)), closes [#159](https://github.com/PharmGKB/PharmCAT/issues/159)
* **reporter:** fix version check for CPIC content ([cb4b762](https://github.com/PharmGKB/PharmCAT/commit/cb4b7627e8010ba463255321915868033d41e007)), closes [#157](https://github.com/PharmGKB/PharmCAT/issues/157)

## [2.8.2](https://github.com/PharmGKB/PharmCAT/compare/v2.8.1...v2.8.2) (2023-10-02)


### Bug Fixes

* **data:** update to PharmVar 6.0.7 ([34e253b](https://github.com/PharmGKB/PharmCAT/commit/34e253b92ee0922162bc2026adee47522a91dac7))
* **namedallelematcher:** fix ArrayIndexOutOfBoundsException bug ([a29420a](https://github.com/PharmGKB/PharmCAT/commit/a29420a43b5ffba505c60fbc42e14fe9de05636e)), closes [#156](https://github.com/PharmGKB/PharmCAT/issues/156)
* **phenotyper:** consistently sort GenePhenotype.diplotypes ([36f1aa5](https://github.com/PharmGKB/PharmCAT/commit/36f1aa509540c9133bef66779c40d23b72e9e659))

## [2.8.1](https://github.com/PharmGKB/PharmCAT/compare/v2.8.0...v2.8.1) (2023-09-29)


### Bug Fixes

* improve error handling in batch mode ([b7c9e66](https://github.com/PharmGKB/PharmCAT/commit/b7c9e660464313b0e42c59e014c357aa25a7455b))
* **data:** add more data validation ([0388049](https://github.com/PharmGKB/PharmCAT/commit/0388049ba06b345f54a4bf2682f872379156d4ce))
* **namedallelematcher:** fix DPYD HapB3 phasing issue ([256ac64](https://github.com/PharmGKB/PharmCAT/commit/256ac64335d47b53c6c7b566b75bc9b0f52a6646)), closes [#155](https://github.com/PharmGKB/PharmCAT/issues/155)

## [2.8.0](https://github.com/PharmGKB/PharmCAT/compare/v2.7.1...v2.8.0) (2023-09-22)


### Features

* **pharmcat:** update DPYD matcher algorithm ([763b950](https://github.com/PharmGKB/PharmCAT/commit/763b950a34f6227eaf759135cbfeec4cd8ea36f5)), closes [#150](https://github.com/PharmGKB/PharmCAT/issues/150)


### Bug Fixes

* **pharmcat:** allow reporter JSON in research mode ([a3b0496](https://github.com/PharmGKB/PharmCAT/commit/a3b04963017e0e99ac3851bd666bdf70513d1262))
* **pharmcat:** avoid NPE ([a86daa6](https://github.com/PharmGKB/PharmCAT/commit/a86daa61e765bb158e7996c0d250bcb1ac7d7304))
* **pharmcat:** avoid NPE ([da11bd6](https://github.com/PharmGKB/PharmCAT/commit/da11bd6556298b109e17543faf00f378cc0b6444))
* **pharmcat:** fix console message ([f2b4086](https://github.com/PharmGKB/PharmCAT/commit/f2b4086b6bfcee20f88874b7c8fea59fc782b8d9))
* **pharmcat:** fix how warnings are handled with outside calls ([634f3b4](https://github.com/PharmGKB/PharmCAT/commit/634f3b4231dd05e1b6cf43eda035ebe7a74afb24)), closes [#154](https://github.com/PharmGKB/PharmCAT/issues/154)
* **phenotyper:** pass through outside calls even if no recommendations are available ([03f7f7c](https://github.com/PharmGKB/PharmCAT/commit/03f7f7cb6f97dab9726f7dc4b180ceb45844fa75)), closes [#154](https://github.com/PharmGKB/PharmCAT/issues/154)
* **reporter:** add DPYD warnings ([2427813](https://github.com/PharmGKB/PharmCAT/commit/24278138e706faa7e749d87debf80b016a4a7683))

## [2.7.1](https://github.com/PharmGKB/PharmCAT/compare/v2.7.0...v2.7.1) (2023-09-09)


### Bug Fixes

* **data:** revert removal of DPYD wobble handling ([9a78c86](https://github.com/PharmGKB/PharmCAT/commit/9a78c86e26bb1dc743a108d14146e7a579e02941))
* **data:** update to PharmVar 6.0.5 ([63c717b](https://github.com/PharmGKB/PharmCAT/commit/63c717bad5c5539b3c9045f047675756e7fd07d6))

## [2.7.0](https://github.com/PharmGKB/PharmCAT/compare/v2.6.0...v2.7.0) (2023-09-02)


### Features

* **pharmcat:** Use PharmGKB for all drug and phenotype annotations ([7390fa1](https://github.com/PharmGKB/PharmCAT/commit/7390fa18fb6f1942266f2278b0852efa1882a46d))
* **reporter:** add new footnote about CPIC/DPWG function ([b986a13](https://github.com/PharmGKB/PharmCAT/commit/b986a138f5ad9266ca3024ee10c9e55a19441b81))


### Bug Fixes

* **data:** fix guideline citation data and update RYR1 phenotypes ([93de1ff](https://github.com/PharmGKB/PharmCAT/commit/93de1ffa7d7f72613dd9de6a48f091744e81dc77))
* **pharmcat:** disable reporter module when using research mode ([d1b1822](https://github.com/PharmGKB/PharmCAT/commit/d1b182225fa4dc033ca2766bf5f7f8513fa1ab00))
* **preprocessor:** capture both error and log from subprocess ([c115ed2](https://github.com/PharmGKB/PharmCAT/commit/c115ed253a489969330afbf56a66b4c2afe5ca39))
* **preprocessor:** install the scikit-allel that has a working toml file ([77e2c43](https://github.com/PharmGKB/PharmCAT/commit/77e2c431f68f2bca3d035b77bff43d012516d15e))

## [2.6.0](https://github.com/PharmGKB/PharmCAT/compare/v2.5.0...v2.6.0) (2023-08-08)


### Features

* **data:** update data to sync with CPIC 1.28 and latest PharmGKB ([4cdfc34](https://github.com/PharmGKB/PharmCAT/commit/4cdfc34d990ebb320a9824831a90d136acfabd7b))


### Bug Fixes

* **data:** update data ([76bfe08](https://github.com/PharmGKB/PharmCAT/commit/76bfe08c85fc1b63d014445e3d7dd883cf7088cf))
* **data:** update F5 version ([43b6637](https://github.com/PharmGKB/PharmCAT/commit/43b6637904c11d2f396c48494ae5089a3ffb36ec))
* **namedallelematcher:** fix DPYD call with ref and partial ([fa09226](https://github.com/PharmGKB/PharmCAT/commit/fa09226f9494d8af680b5cad442d174dfbef9b5e))
* **namedallelematcher:** improve DPYD calling ([57b2bd8](https://github.com/PharmGKB/PharmCAT/commit/57b2bd81d38f6ac871ad9a4e410df700c4b11def))
* **namedallelematcher:** improve DPYD calling with HapB3 wobble ([3528815](https://github.com/PharmGKB/PharmCAT/commit/3528815318faf9d49370cae97c96c50db62e90d6))
* **namedallelematcher:** support wobble in DPYD HapB3 ([65a1309](https://github.com/PharmGKB/PharmCAT/commit/65a13098f46a43d422058a34d48a2667bed184c3))
* **namedallelematcher:** update HTML output to use Bootstrap 5 ([3774cca](https://github.com/PharmGKB/PharmCAT/commit/3774cca4108268e71d4c811546db3b46fedabd82))
* **pharmcat:** fix NPE when called without -o flag ([a8d5f51](https://github.com/PharmGKB/PharmCAT/commit/a8d5f5100faa5560a8ffefbf56c67e9a27f1029f))
* **preprocessor:** pipe stdout to a user's running interface to show warning messages from a subprocess ([36c151b](https://github.com/PharmGKB/PharmCAT/commit/36c151b568392942065d8b41d3be966df399c503))
* **reporter:** display allele function in section 3 ([5f48358](https://github.com/PharmGKB/PharmCAT/commit/5f48358447780ae24cb920fd2698655aff6a2b0f))
* **reporter:** don't emit empty ids ([c7f3d48](https://github.com/PharmGKB/PharmCAT/commit/c7f3d480c3125fb344569ec151b3888f49166aa0))
* **reporter:** fix browser upgrade CSS ([6ba5bf1](https://github.com/PharmGKB/PharmCAT/commit/6ba5bf1d5341d64250fe6d1212e39f5fec1f6dca))
* **reporter:** fix links to uncallable genes in section 3 ([f1d4da2](https://github.com/PharmGKB/PharmCAT/commit/f1d4da2c46fcb0aa6a87643e3a0d976f215992e3))
* **reporter:** identify homozygous DPYD haplotypes ([0f10fe3](https://github.com/PharmGKB/PharmCAT/commit/0f10fe386ae2cccc6bd64b770dcb896a9c3f9c24))
* **reporter:** improve allele match data for CYP2C19 *1 ([2bc80cd](https://github.com/PharmGKB/PharmCAT/commit/2bc80cdd31d0bd8818c542037a578769a71e2583))
* **reporter:** improve message for uncalled because of no data ([8211015](https://github.com/PharmGKB/PharmCAT/commit/82110153b980b2eb25b6b3bdc9131fdc99f19a35))
* **reporter:** improve text for uncallable genes ([7e794b3](https://github.com/PharmGKB/PharmCAT/commit/7e794b380bfb9c77b67d289abeb204526cd31698))
* **reporter:** list function in section 3 ([e51483c](https://github.com/PharmGKB/PharmCAT/commit/e51483ca14dfabdfa7257faa0292c199699f00c3))
* **reporter:** show function in section 3 even when there is a no call ([4e2ea20](https://github.com/PharmGKB/PharmCAT/commit/4e2ea20fd4165918a538707407a2f0f7a44abaf8))
* **reporter:** show genes in section III if uncallable ([7da2075](https://github.com/PharmGKB/PharmCAT/commit/7da2075732beecb856430dc06ac45356cbc47e1b))

## [2.5.0](https://github.com/PharmGKB/PharmCAT/compare/v2.4.0...v2.5.0) (2023-06-12)


### Features

* **namedallelematcher:** support treating undocumented variation as reference ([f5a5baa](https://github.com/PharmGKB/PharmCAT/commit/f5a5baa5b984f580e1763c2028649148bd268e00))


### Bug Fixes

* fix loading of path for test file writing ([beb03c3](https://github.com/PharmGKB/PharmCAT/commit/beb03c337944458fe86dee7e3d11d68a39a48943))
* **data:** update CPIC data ([ab0fb0f](https://github.com/PharmGKB/PharmCAT/commit/ab0fb0f2e56c2d4133058ea9360bea5366e91890))
* **namedallelematcher:** adhere to AD format rules ([561012a](https://github.com/PharmGKB/PharmCAT/commit/561012acec6722c47414e04e6f93f5d87adf8726)), closes [#139](https://github.com/PharmGKB/PharmCAT/issues/139)
* **reporter:** display consistent genotype in report, fix SLCO1B1 bug ([08f65ea](https://github.com/PharmGKB/PharmCAT/commit/08f65ea2c345331c1167cf9625227be2850b449d))
* **reporter:** fix DPWG version check ([6c2d828](https://github.com/PharmGKB/PharmCAT/commit/6c2d8288bdad8eeab6842e204c060589bbcd362d))
* **reporter:** fix JSON property name ([e7036c8](https://github.com/PharmGKB/PharmCAT/commit/e7036c848eaab8035b39a9681b3c3932ca9803aa))
* **reporter:** improve styling when there are multiple drug recommendations ([20b20f4](https://github.com/PharmGKB/PharmCAT/commit/20b20f435528e012c97c7ff0dd98210a85f7b05e))

## [2.4.0](https://github.com/PharmGKB/PharmCAT/compare/v2.3.0...v2.4.0) (2023-05-01)


### ⚠ BREAKING CHANGES

* update to yarn 3.5

### Features

* update data for CPIC and PharmVar content ([8fe4173](https://github.com/PharmGKB/PharmCAT/commit/8fe41737457e8f755f844b4adff31baceb1b2b14))
* update to latest CPIC/PharmVar data for CYP2D6 and SSRI/SNRIs ([ebaea81](https://github.com/PharmGKB/PharmCAT/commit/ebaea81c5d807866219e349f7a005130fdbe1ce5))


### Bug Fixes

* improve how DefinitionReader is instantiated ([8596063](https://github.com/PharmGKB/PharmCAT/commit/85960633a6887bb74d3a6d460a57b45f9fb3aee0))
* standardize output order ([df675c2](https://github.com/PharmGKB/PharmCAT/commit/df675c21b7520ae9d5967010269f41f243129e29))
* **data:** update messages ([2a371dc](https://github.com/PharmGKB/PharmCAT/commit/2a371dc41322b03fe6886b795caf53e1218e10bd))
* **data:** update messages ([1b4dd8b](https://github.com/PharmGKB/PharmCAT/commit/1b4dd8bae56c3a977effbc2c1aabdfb0e980e442))
* **namedallelematcher:** do not expose novelAllele ([cc93ed0](https://github.com/PharmGKB/PharmCAT/commit/cc93ed0b21ef81dafc58310b755d8818676db8cb))
* **namedallelematcher:** fix null in toString ([cbf827a](https://github.com/PharmGKB/PharmCAT/commit/cbf827a4a4e115bccaf4dc5eb21842fd607a9d99))
* **namedallelematcher:** track novel alleles for sample and fix erroneous novel allele error when novel allele occurs in different sample ([b079c54](https://github.com/PharmGKB/PharmCAT/commit/b079c542f0a0a394f7293725ec2a8d98cbbf48b8))
* **report:** remove max-width from print style ([d5909e0](https://github.com/PharmGKB/PharmCAT/commit/d5909e0d5df7b03749d79c481eda816dbe8e9919)), closes [#129](https://github.com/PharmGKB/PharmCAT/issues/129)
* **reporter:** handle missing gene report ([f9d7dba](https://github.com/PharmGKB/PharmCAT/commit/f9d7dba7892b0c44a8d521fa14b5d732ca6d8a7e))
* **reporter:** show gene in section 3 if uncallable but has data ([264cb0e](https://github.com/PharmGKB/PharmCAT/commit/264cb0eae53b175a62c7b5207648af27d595c7f8))


### Build System

* update to yarn 3.5 ([7d0eeb1](https://github.com/PharmGKB/PharmCAT/commit/7d0eeb1f1400da38fab8b7b0b1560adbf45b4658))

## [2.3.0](https://github.com/PharmGKB/PharmCAT/compare/v2.2.3...v2.3.0) (2023-03-14)


### Features

* add link to license on the homepage ([e71e7c1](https://github.com/PharmGKB/PharmCAT/commit/e71e7c1c379e7f618882a4843b7f88f07f3b0115))
* improve styling for print media ([c94ca95](https://github.com/PharmGKB/PharmCAT/commit/c94ca95e881a8fa936694d30f02bb1ad9e8e3180))


### Bug Fixes

* bad link on homepage ([03c788d](https://github.com/PharmGKB/PharmCAT/commit/03c788d9ae9a909fad9e167c956db3bc828f7f0c)), closes [#127](https://github.com/PharmGKB/PharmCAT/issues/127)
* update data ([c526cfc](https://github.com/PharmGKB/PharmCAT/commit/c526cfcd781c439efc54ef4d61006c551d651b9f))
* update data ([209032d](https://github.com/PharmGKB/PharmCAT/commit/209032d03b03d57e4cdabc531bd33debb0c26b6d))
* **namedallelematcher:** improve error message when invalid GT allele value is provided in VCF ([4af5de0](https://github.com/PharmGKB/PharmCAT/commit/4af5de00f926b980999e2ea6181294d122be9f81))
* **pipeline:** support specifying max java heap size for pipeline ([bdea14f](https://github.com/PharmGKB/PharmCAT/commit/bdea14f5c3d9f1c86fa655e9fd095f0f25ebe30b))
* **reporter:** wrap long calls nicely ([9bc9cfb](https://github.com/PharmGKB/PharmCAT/commit/9bc9cfbcbadbfc64b1a01e246f485774c96646b9)), closes [#130](https://github.com/PharmGKB/PharmCAT/issues/130)

## [2.2.3](https://github.com/PharmGKB/PharmCAT/compare/v2.2.2...v2.2.3) (2023-02-06)


### Bug Fixes

* fix permissions on executables ([3385186](https://github.com/PharmGKB/PharmCAT/commit/3385186c3d8e70cf7cd0c34511878e3650b96da0))
* update to latest CPIC/PharmVar data ([905306f](https://github.com/PharmGKB/PharmCAT/commit/905306ffa2b5c932f8ddedc238dd1593ed65c51f))
* **reporter:** fix case for allele function names ([90922b5](https://github.com/PharmGKB/PharmCAT/commit/90922b59d444da7257eb3f97f6c0a8ecd1d8bcf8))

## [2.2.2](https://github.com/PharmGKB/PharmCAT/compare/v2.2.1...v2.2.2) (2023-02-01)


### Bug Fixes

* update to latest data from PharmVar and CPIC ([6ac3aeb](https://github.com/PharmGKB/PharmCAT/commit/6ac3aebb4719379f77f10137d0933e3127959814)), closes [#126](https://github.com/PharmGKB/PharmCAT/issues/126)
* **docker:** sanitize ownership of reference FASTA files ([2599f4e](https://github.com/PharmGKB/PharmCAT/commit/2599f4e3e000474c96108d151a825f7360ed310c))
* **pharmcat:** add pharmcat_pipeline script ([5c2fd85](https://github.com/PharmGKB/PharmCAT/commit/5c2fd851bde215c6e527491ce0331a61e03371e7))
* **pharmcat:** fix max -cp recommendation ([76295d0](https://github.com/PharmGKB/PharmCAT/commit/76295d08dc8e5fbc95f4cbe5e25b1595e0418f43))

## [2.2.1](https://github.com/PharmGKB/PharmCAT/compare/v2.2.0...v2.2.1) (2023-01-21)


### Bug Fixes

* **preprocessor:** strip AN/AC info out of single sample output ([8554572](https://github.com/PharmGKB/PharmCAT/commit/8554572e67ae7402cbb833a3e05242fc3bbc954d))

## [2.2.0](https://github.com/PharmGKB/PharmCAT/compare/v2.1.2...v2.2.0) (2023-01-20)


### Features

* **namedallelematcher:** support multisample VCF ([4790ac0](https://github.com/PharmGKB/PharmCAT/commit/4790ac02fb4b3319171f51cb569d71982d394d60))
* **pharmcat:** add concurrency support to pharmcat ([eba4c4c](https://github.com/PharmGKB/PharmCAT/commit/eba4c4cae7c0e99762fe52724aeffacd90864682))


### Bug Fixes

* **data:** update data ([77e5430](https://github.com/PharmGKB/PharmCAT/commit/77e54303bc6571070e44354e8374a1ece9802466))
* **namedAlleleMatcher:** fix handling of homozygous combinations ([e13b587](https://github.com/PharmGKB/PharmCAT/commit/e13b58704cd436b2f3dad6681b7af7ea1d5967c7))
* **namedallelematcher:** fix handling of homozygous result for DPYD ([ea127ea](https://github.com/PharmGKB/PharmCAT/commit/ea127eaedb2ad59fbdd12939a61547929ea3fd9b))
* **namedallelematcher:** fix parsing of AD field ([b6b0d73](https://github.com/PharmGKB/PharmCAT/commit/b6b0d7346870ed26cfa0a8e65e82d0ae8dc7e36e))
* **namedallelematcher:** fix parsing of AD field ([04326ef](https://github.com/PharmGKB/PharmCAT/commit/04326ef994d0064e91e3a880f9284eaa97ac17b9)), closes [#118](https://github.com/PharmGKB/PharmCAT/issues/118)
* **namedallelematcher:** fix preprocessed spelling ([e7c4a65](https://github.com/PharmGKB/PharmCAT/commit/e7c4a65e05b82552b5047a29c42f554824390249))
* **namedallelematcher:** haplotype name comparison ([8e1f547](https://github.com/PharmGKB/PharmCAT/commit/8e1f5479301486e5212766f53347b45e844e7150)), closes [#119](https://github.com/PharmGKB/PharmCAT/issues/119)
* **namedallelematcher:** read VCF into memory if possible ([e1bc46d](https://github.com/PharmGKB/PharmCAT/commit/e1bc46d256c2a4baabd7e56cfcce8b98b8305b66))
* **namedallelematcher:** support vcf.gz and vcf.bgz ([1710217](https://github.com/PharmGKB/PharmCAT/commit/1710217a1c2551f71de4254d76444c7aa6489ce5))
* **pharmcat:** fix -def flag parsing ([3407310](https://github.com/PharmGKB/PharmCAT/commit/34073104610671ad10875292a27237b4cdffef71))
* **pharmcat:** fix BatchPharmCAT problems ([7c045ff](https://github.com/PharmGKB/PharmCAT/commit/7c045ff0607b2f96e4f371b7693166f633b701d1))
* **pharmcat:** improve pipeline, carry sample id throughout ([1a3b954](https://github.com/PharmGKB/PharmCAT/commit/1a3b9548bb0eb4bd19859e966d473006446cd342))
* **pharmcat:** support specifying sample; sort data model to keep tests stable ([05007e3](https://github.com/PharmGKB/PharmCAT/commit/05007e3ac76e7533c7ca8cf33bf077ff992e2c26))
* **preprocessor:** add java check, support tool paths via environment variables ([a762127](https://github.com/PharmGKB/PharmCAT/commit/a76212754de711550a523ea703663dd3d06ace0c))
* **preprocessor:** change the way to determine positions without any called genotype ([ea43698](https://github.com/PharmGKB/PharmCAT/commit/ea43698aa89090d82d104f80a7e8dccd28e44a51))
* **preprocessor:** fix bugs introduced by DPYD change ([c748725](https://github.com/PharmGKB/PharmCAT/commit/c7487259f75e229798ba82b47d737d77b4a9c283))
* **preprocessor:** fix output filename handling ([ec26f66](https://github.com/PharmGKB/PharmCAT/commit/ec26f66501312969faaed6bbf0babd2d908c4029))
* **preprocessor:** move gvcf check to beginning ([c9d3837](https://github.com/PharmGKB/PharmCAT/commit/c9d38374f34429d977f9421153c37b2935e5dfe8))
* **preprocessor:** output single VCF by default ([2bf0c74](https://github.com/PharmGKB/PharmCAT/commit/2bf0c7415ec5629cd14e1c8940f2c4c5ddbc20ed))
* **preprocessor:** print out the tool version ([9edc5ef](https://github.com/PharmGKB/PharmCAT/commit/9edc5ef26578f1a8f0672f1a2d62a6fce5186bda))
* **reporter:** rename "generatedOn" property to "timestamp" ([c241d01](https://github.com/PharmGKB/PharmCAT/commit/c241d0125354b224cdb807a3744bc5189ae031f3))

## [2.1.2](https://github.com/PharmGKB/PharmCAT/compare/v2.1.1...v2.1.2) (2022-11-11)


### Bug Fixes

* **namedAlleleMatcher:** add support for PCATxINDEL from preprocessor ([bc04d5d](https://github.com/PharmGKB/PharmCAT/commit/bc04d5d41c8a5a12bffd9152b7943ee43a2a769a))
* **preprocessor:** bug in determining indel ([35e25cc](https://github.com/PharmGKB/PharmCAT/commit/35e25cc98efd182ae62d525fd0da7bd39860d477))
* **preprocessor:** fix a bug related to missing genotypes ([ff676c3](https://github.com/PharmGKB/PharmCAT/commit/ff676c30b73ae21e5969599fdd8907a6b80e910f))
* **preprocessor:** fix a bug with the unspecified ALT ([88b94d2](https://github.com/PharmGKB/PharmCAT/commit/88b94d289d403a6ff14091d9cf88bb4f9e2f6257))
* **preprocessor:** fix error for indels ([0bc0b74](https://github.com/PharmGKB/PharmCAT/commit/0bc0b7432f83c1db529e15908cb262ba509e3e4c))
* **preprocessor:** fix pharmcat positions file lookup logic ([9f2a77a](https://github.com/PharmGKB/PharmCAT/commit/9f2a77a7b2de79bf7ef03d1517cd0a42d69e116e))

## [2.1.1](https://github.com/PharmGKB/PharmCAT/compare/v2.1.0...v2.1.1) (2022-11-09)


### Bug Fixes

* **preprocessor:** add version info and improve docs ([8d1d785](https://github.com/PharmGKB/PharmCAT/commit/8d1d785351b4d7e9865eeec6e812a210ff9f62d8))
* **preprocessor:** fix automatic version update ([c3d5940](https://github.com/PharmGKB/PharmCAT/commit/c3d5940c2b778849c3859a26006022f5d28b4fff))

## [2.1.0](https://github.com/PharmGKB/PharmCAT/compare/v2.0.3...v2.1.0) (2022-11-09)


### Features

* **preprocessor:** major refactoring and support for concurrent mode ([ce5c3e5](https://github.com/PharmGKB/PharmCAT/commit/ce5c3e5b9a5d58d8a3a677f0464246cd014c9201))


### Bug Fixes

* bug in reference genome sequence ([6fc5550](https://github.com/PharmGKB/PharmCAT/commit/6fc555015cc997af018656266787133d5db9724e))
* **data:** update to CPIC 1.21.3 ([61143e6](https://github.com/PharmGKB/PharmCAT/commit/61143e68200da9a904791cff96f6a6f67b0026b6))
* **preprocessor:** fix a bug for alt=<*> ([c3b96c3](https://github.com/PharmGKB/PharmCAT/commit/c3b96c3969fd208f0a26eca32dc4fe9ad39558c6))
* **preprocessor:** process regions concurrently ([7c4f581](https://github.com/PharmGKB/PharmCAT/commit/7c4f581c1381962477c9753ee08a82d30bb1b47a))


### Performance Improvements

* **preprocessor:** add support for concurrent mode ([71855f2](https://github.com/PharmGKB/PharmCAT/commit/71855f2117599673a4678cfd91a5a93fc873e7d1))

## [2.0.3](https://github.com/PharmGKB/PharmCAT/compare/v2.0.2...v2.0.3) (2022-10-27)


### Bug Fixes

* **pharmcat:** fix handling of relative files ([fa8e76f](https://github.com/PharmGKB/PharmCAT/commit/fa8e76fc51f83e201e7565042081a140e7c199fe))
* **preprocessor:** fix preprocessor bug ([d542e52](https://github.com/PharmGKB/PharmCAT/commit/d542e520b283a19e3eac5ae74b4665dd8fd6ad12))

## [2.0.2](https://github.com/PharmGKB/PharmCAT/compare/v2.0.1...v2.0.2) (2022-10-26)


### Bug Fixes

* errors in the names of sex and mitochondria chromosomes in the reference genome sequence ([b5e1412](https://github.com/PharmGKB/PharmCAT/commit/b5e1412c5cab74eb5b046b195302c914d7259293)), closes [#116](https://github.com/PharmGKB/PharmCAT/issues/116)
* fix support for haploid genes (chrX, chrY, chrM) ([dc8a43d](https://github.com/PharmGKB/PharmCAT/commit/dc8a43df69cda82957badbf88cd62297c4c1fab7))
* **data:** fix problem with chrX data preprocessor relies on ([2d71671](https://github.com/PharmGKB/PharmCAT/commit/2d716716dd7b8e6636fe43d619b13cb78f3a646e))
* **data:** update to latest CPIC data ([4393cb3](https://github.com/PharmGKB/PharmCAT/commit/4393cb33dc58b8361b24274c2dbc32146dd084f5))

## [2.0.1](https://github.com/PharmGKB/PharmCAT/compare/v2.0.0...v2.0.1) (2022-10-25)


### Bug Fixes

* **preprocessor:** declare and concat list variables ([d0899f0](https://github.com/PharmGKB/PharmCAT/commit/d0899f01e2139b91b6deb081adeaceee235ca5a6)), closes [#115](https://github.com/PharmGKB/PharmCAT/issues/115)
* **preprocessor:** fix a bug in bcftools and bgzip version checks that interrupts the VCF preprocessor ([34f3f54](https://github.com/PharmGKB/PharmCAT/commit/34f3f54d645081042e32cd86b20d4abe127f7dc8)), closes [#116](https://github.com/PharmGKB/PharmCAT/issues/116)
* **preprocessor:** fix a bug that interrupts the VCF preprocessor ([0ed6e22](https://github.com/PharmGKB/PharmCAT/commit/0ed6e2233859a1f33950bbe1694c4ca4be6c4c46)), closes [#115](https://github.com/PharmGKB/PharmCAT/issues/115)

## [2.0.0](https://github.com/PharmGKB/PharmCAT/compare/v1.8.0...v2.0.0) (2022-10-20)


### ⚠ BREAKING CHANGES

* **reporter:** this introduces changes to the phenotyper and reporter JSON
files that are not backwards-compatible.
* org.pharmgkb.pharmcat.PharmCAT is now the only command line tool.
All arguments have been changed except for -vcf.

This change will only impact you if you specify anything other than -vcf or use
anything other than the main PharmCAT program.
* The command line arguments to modules have been standardized to match the main PharmCAT program.
  * NamedAlleleMatcher: -d renamed -na, -a renamed -ar
  * Phenotyper: -o renamed -a
  * Reporter: -o renamed -f
* **namedallelematcher:** The NamedAlleleMatcher will always assume reference.  This is no longer a runtime option.
* PharmCAT now requires Java 17

### Features

* gate CYP2D6 calling behind research flag ([12e4fdf](https://github.com/PharmGKB/PharmCAT/commit/12e4fdfba822d5953b0c4ea33daa1a364f760522))
* require Java 17 ([5ca456e](https://github.com/PharmGKB/PharmCAT/commit/5ca456e8ea92eded6b45c01570e017cd70a5f566))
* update PharmCAT pipeline tool arguments ([08e6e22](https://github.com/PharmGKB/PharmCAT/commit/08e6e22eb4278d0b8bcd65da238bda3596538ea0))
* **data:** add definition for CYP2D6 ([f0ab563](https://github.com/PharmGKB/PharmCAT/commit/f0ab56321a32afa4b745d051d52d6a4adcccafdd))
* **datamanager:** add DPWG allele defintions including CYP3A4 ([1842f32](https://github.com/PharmGKB/PharmCAT/commit/1842f3207de8ce92dc30e9fe85ea8fc1e08b82d7))
* **datamanager:** add DPWG file processing to the DataManager ([58f7691](https://github.com/PharmGKB/PharmCAT/commit/58f769146f882d60bd5476e29076b47f45bb8204))
* **namedallelematcher:** add support for combinations and partials, update version to 2.0 ([d3a6065](https://github.com/PharmGKB/PharmCAT/commit/d3a6065bf29c09ad6894c635d1080b64f518fd31))
* **namedallelematcher:** removing assumeReference ([864fea0](https://github.com/PharmGKB/PharmCAT/commit/864fea097de83e684357f2abb15ae98b7300bda8))
* **phenotyper:** add activity data to diplotypes and haplotypes ([ec12285](https://github.com/PharmGKB/PharmCAT/commit/ec122850120b4a2bd783aa17088cec253b1ced90))
* **phenotyper:** make sample VCF data take precedence over outside calls ([dc8f585](https://github.com/PharmGKB/PharmCAT/commit/dc8f585cd2d444f1725d931352e30d2f32eabee4))
* **phenotyper:** support outside activity score ([52cfe3c](https://github.com/PharmGKB/PharmCAT/commit/52cfe3ca42260e43b1b504e13cdba1e4db9fdfa4))
* **phenotyper:** switch to preferring outside call data over VCF data and update wording accordingly ([a3451aa](https://github.com/PharmGKB/PharmCAT/commit/a3451aae557efab2f0cf5582436c21854036299b))
* **reporter:** add back comments and activity score to report ([1f55f7c](https://github.com/PharmGKB/PharmCAT/commit/1f55f7cd4c8948fb3c560902e11acba987032505))
* **reporter:** add DPWG drug data ([7a32257](https://github.com/PharmGKB/PharmCAT/commit/7a32257241e16aedba9548346702cf586acfb425))
* **reporter:** add DPYD-specific allele matcher for use in unphased data ([2678251](https://github.com/PharmGKB/PharmCAT/commit/2678251cacb2071cda3eb747026233853ebafbe8))
* **reporter:** add messages for combo and cyp2d6 modes ([9604086](https://github.com/PharmGKB/PharmCAT/commit/960408645e541d53abc93ad82191cee821d208e7))
* **reporter:** Add new Genotype, GuidelineReport, and AnnotationGroup models for organizing multiple results ([c468223](https://github.com/PharmGKB/PharmCAT/commit/c468223b73e3b8dc48025301791118f423475439))
* **reporter:** add source labels to gene and drug summary ([c6bda14](https://github.com/PharmGKB/PharmCAT/commit/c6bda14836282bc1fa82ea92ee3ec9b79838d639))
* **reporter:** enable DrugCollection and GeneDrugSummary to run independently ([17c773d](https://github.com/PharmGKB/PharmCAT/commit/17c773dddd6e62c5558bedd60374886eedb4942b))
* **reporter:** fix genotype matching for DPWG guidelines ([e66d8c0](https://github.com/PharmGKB/PharmCAT/commit/e66d8c0b8559a9be03ba9c50719ba63c21affca0))
* **reporter:** group genotype, function, and phenotype in gene summary of HTML report ([d885cf0](https://github.com/PharmGKB/PharmCAT/commit/d885cf0445d16634913cc083cb9a607437f5d93f))
* **reporter:** label inferred diplotype calls in the final report data ([32c48dc](https://github.com/PharmGKB/PharmCAT/commit/32c48dc7490a66882d0ef7a3ab8ba298850bed5d))
* **reporter:** make least-function allele calling work for DPYD in the reporter ([5820d80](https://github.com/PharmGKB/PharmCAT/commit/5820d808bbfc014862f903dab0606e1d887fd50a))
* **reporter:** more work to change to GuidelineReport for data tracking ([bf2fc7e](https://github.com/PharmGKB/PharmCAT/commit/bf2fc7ed96410100da3d1ddc87003017e8b48a0e))
* **reporter:** move population into drug table row header and fix row styling ([2f610d2](https://github.com/PharmGKB/PharmCAT/commit/2f610d23a4702d16ad637261b7b2215f9f73d881))
* **reporter:** rename "CPIC" recommendations to "Prescribing" recommendations ([9d10303](https://github.com/PharmGKB/PharmCAT/commit/9d10303ea9e60afbf49f0648a548c0117a85e061))
* **reporter:** reorder implications and recommendations columns in final report ([e7542a1](https://github.com/PharmGKB/PharmCAT/commit/e7542a1cf41878863df2efa19fa2202dc56cd627))
* **reporter:** separate CPIC/DPWG reports ([87e4bd2](https://github.com/PharmGKB/PharmCAT/commit/87e4bd21009e4894520e322ab85a4ebad5f53c54))
* **reporter:** split Reporter class into separate AbstractFormat classes and retain ReportContext ([d2ebd53](https://github.com/PharmGKB/PharmCAT/commit/d2ebd53264f3e99cb4fbfb68173d007659fe9eb2))
* **reporter:** switch DPWG guideline data source to individual files ([bd0e180](https://github.com/PharmGKB/PharmCAT/commit/bd0e1800947f28e7b15d6dbbe8e1be6c9243e578))
* **reporter:** update formatting of the HTML report ([a65913c](https://github.com/PharmGKB/PharmCAT/commit/a65913c7fa67a1116f8afce0891e937b70f2b56c))
* **website:** add a make task to publish a prototype of the PharmCAT website ([5b8383b](https://github.com/PharmGKB/PharmCAT/commit/5b8383b685f484a7c031a0303661cbbaee3e2d78))
* **reporter:** improve function display for DPYD combo genotypes ([b8a29ef](https://github.com/PharmGKB/PharmCAT/commit/b8a29ef0145b4c861701bfd589d6592224a4a32b))


### Bug Fixes

* normalize cli args ([f754e50](https://github.com/PharmGKB/PharmCAT/commit/f754e50c6ca1242371f1172f66879ea77a6464a8))
* revert to GRCh38.p13 ([0b0eeda](https://github.com/PharmGKB/PharmCAT/commit/0b0eeda45074abd522849907244070239de3128c))
* standardize source/version info throughout pipeline ([679842c](https://github.com/PharmGKB/PharmCAT/commit/679842cc3021f0132674d2e504554201b5e776cc))
* update DataManager ([61fb5d1](https://github.com/PharmGKB/PharmCAT/commit/61fb5d1dfa3d0733634766b710659bab5c7e6304))
* update to GRCh38.p14 ([8fe2728](https://github.com/PharmGKB/PharmCAT/commit/8fe27287af5ad4ca7be3e5f4664d418118de78fa))
* update to most recent PharmVar/CPIC/DPWG ([c689178](https://github.com/PharmGKB/PharmCAT/commit/c68917815b4a5e06ca0adf1b74b874207a84cc0a))
* make command line flags more consistent ([28b31fb](https://github.com/PharmGKB/PharmCAT/commit/28b31fbe78deb07b00f6d952f179578f6ddc8d9a))
* **data:** add support for G6PD ([5790da6](https://github.com/PharmGKB/PharmCAT/commit/5790da62b729c2d0f478a969a02ea07d1cccaf0b))
* **data:** remove "total activity score" in DPWG data ([db3273e](https://github.com/PharmGKB/PharmCAT/commit/db3273ec04f5b0751cd4ad67bae3418941fda5ab))
* **data:** update data ([3898afe](https://github.com/PharmGKB/PharmCAT/commit/3898afe0d13dd1d12bb0be69f0280ad79c216ca6))
* **data:** update data ([721f77a](https://github.com/PharmGKB/PharmCAT/commit/721f77a07bf00f07ee30b021a453148360d8bef4))
* **data:** update data ([fe8c19a](https://github.com/PharmGKB/PharmCAT/commit/fe8c19a2fa82a7fdaed6cfea2142bfcdeedbeb70))
* **data:** update data ([d378fc8](https://github.com/PharmGKB/PharmCAT/commit/d378fc85bd78fc30e0e11c81084951a61e02c011))
* **data:** update data to latest CPIC/DPWG ([5f2cedb](https://github.com/PharmGKB/PharmCAT/commit/5f2cedbf08b0bb0e379bcbf410cdc136628fe7fd))
* **data:** update messages ([3780c7d](https://github.com/PharmGKB/PharmCAT/commit/3780c7d1d1cdbb3f1508d456b35e1875790c93a2))
* **datamanager:** fix exemptions for ignored positions/alleles ([ee9a41e](https://github.com/PharmGKB/PharmCAT/commit/ee9a41ecca99b08660fc83a18cc4e5f994e88a19))
* **datamanager:** fix logic and use proper exit code on error ([1ca34a0](https://github.com/PharmGKB/PharmCAT/commit/1ca34a0bbaff7fdaf686f20d9f3ead33a10e2a56))
* **datamanager:** handle 429 error codes ([28b4422](https://github.com/PharmGKB/PharmCAT/commit/28b4422ebca1baff6caa9fde5007c1722f7290b6))
* **datamanager:** improve -sdl handling ([c72a2e8](https://github.com/PharmGKB/PharmCAT/commit/c72a2e82efffe88b1100d56c2d205f25a2689734))
* **datamanager:** support HGVS-based positions in exemptions ([584d552](https://github.com/PharmGKB/PharmCAT/commit/584d552118ade52e84e270e70dfaa0ce12d2b2fa))
* **docker:** update samtools versions ([2486315](https://github.com/PharmGKB/PharmCAT/commit/24863159bb85e12579f80046ab8dee12a816fbf6))
* **namedallelematcher:** add effectively phased property ([397ad5c](https://github.com/PharmGKB/PharmCAT/commit/397ad5cf82f1882d7a45dc3ab5fd0bb5db2b6c66))
* **namedallelematcher:** add runtime metadata ([d988a77](https://github.com/PharmGKB/PharmCAT/commit/d988a77914be00a9434fec9b0d7f55938f06b53f))
* **namedallelematcher:** add wobble support ([1200b1f](https://github.com/PharmGKB/PharmCAT/commit/1200b1fb5bcd569b321b4e36a514e627ffdf38c0))
* **namedallelematcher:** call haplotypes for DPYD if no exact (effectively phased) diplotype match ([492be33](https://github.com/PharmGKB/PharmCAT/commit/492be3367247c99b4d3dd002a6459d8a6eda21cb))
* **namedallelematcher:** fix handling of DPYD Reference ([aa4d095](https://github.com/PharmGKB/PharmCAT/commit/aa4d095af289da4a5f1e0cb88e57bd77669eeaae))
* **namedallelematcher:** fix hetero pair matching in combinations ([9d45ee3](https://github.com/PharmGKB/PharmCAT/commit/9d45ee38363eacf1278acab8b44a3a503fe311e5))
* **namedallelematcher:** improve combination scoring and support off-reference partials ([439b6ba](https://github.com/PharmGKB/PharmCAT/commit/439b6bab1fbd286f2609683d4f115a63ba62ea58))
* **namedallelematcher:** improve DPYD support ([1adbec1](https://github.com/PharmGKB/PharmCAT/commit/1adbec127e7c488d5854065ced2b5abee20ad361))
* **namedallelematcher:** improve matcher html output ([93e75df](https://github.com/PharmGKB/PharmCAT/commit/93e75df9e16ed82bc6c44f2b36caf837d6006f82))
* **namedallelematcher:** properly report phased status ([b774b59](https://github.com/PharmGKB/PharmCAT/commit/b774b59cce1d026b0627f50f88b76cc3ed37e3dc))
* **namedallelematcher:** sort variants by VCF position ([7a2980a](https://github.com/PharmGKB/PharmCAT/commit/7a2980a6a3c42e2f211b70be20a8b2c3731e5f3f))
* **namedallelematcher:** update combination name style ([730af57](https://github.com/PharmGKB/PharmCAT/commit/730af5725c7004ee75c3fd4f72d2137d602ed71a))
* **namedallelematcher:** when assuming reference, use reference for position if reference named allele has wobble ([aca267e](https://github.com/PharmGKB/PharmCAT/commit/aca267efd961876cc5b3bc02debccdc18fc08a2a))
* **pharmcat:** add batch cli ([f25fe58](https://github.com/PharmGKB/PharmCAT/commit/f25fe587b9c7539400d12dfd0996d5d5ccfacfcf))
* **pharmcat:** cannot call System.exit or tests will fail ([789da71](https://github.com/PharmGKB/PharmCAT/commit/789da71a5f5a1a3538cf0524a039fd6edcfdc2d5))
* **pharmcat:** default to compact mode ([6989be0](https://github.com/PharmGKB/PharmCAT/commit/6989be02de6867637dc2acbc85ab2e1bee79f3c0))
* **pharmcat:** exit with proper exit codes ([6b306b6](https://github.com/PharmGKB/PharmCAT/commit/6b306b61401ec4b804346ea672bf1b6522006f72))
* **pharmcat:** fix missing variant warnings when using serialized matcher data ([ed34e28](https://github.com/PharmGKB/PharmCAT/commit/ed34e28b6dfbfd41ae28573c61e4b3f857b6d763))
* **pharmcat:** maintain consistent base output file name ([8ba7b6e](https://github.com/PharmGKB/PharmCAT/commit/8ba7b6e82466b2ab35288ab0bcbc022d6fb98787)), closes [#113](https://github.com/PharmGKB/PharmCAT/issues/113)
* **phenotyper:** add support for calling CYP2D6 ([aaa598a](https://github.com/PharmGKB/PharmCAT/commit/aaa598a4bc4d3f949b459c37754e8cd36814984d))
* **phenotyper:** allow diplotype and phenotype data in the outside calls ([03dd114](https://github.com/PharmGKB/PharmCAT/commit/03dd1146f6daea9e0f5436bf761769c537e81a65))
* **phenotyper:** avoid NPE when applying messages ([fd07844](https://github.com/PharmGKB/PharmCAT/commit/fd07844f25dc26c727c2265883be9d1c547b723d))
* **phenotyper:** fix handling replacement of outside call for GeneReports in Phenotyper ([c62d888](https://github.com/PharmGKB/PharmCAT/commit/c62d888cee391b5ea0693ab3c7afb31d7377e5cb))
* **phenotyper:** fix how DPYD phased alleles are reduced for lookup in Reporter ([25493be](https://github.com/PharmGKB/PharmCAT/commit/25493be0c6a0a2b090959d669105821c88f2edd5))
* **phenotyper:** reorganize OutsideCall into phenotype package ([5cbeb14](https://github.com/PharmGKB/PharmCAT/commit/5cbeb14ca2589064533c73f5e7f1587c1eebf9bb))
* **phenotyper:** standardize GSON usage ([79407e1](https://github.com/PharmGKB/PharmCAT/commit/79407e10edd098d14c2bc46fb592ba0d215fa7d6))
* **preprocessor:** correct phasing status for multiallelic positions ([59fb6b8](https://github.com/PharmGKB/PharmCAT/commit/59fb6b8f9cedf7bd86d0febb86baa8afc8551141)), closes [#102](https://github.com/PharmGKB/PharmCAT/issues/102)
* **preprocessor:** correct phasing status for multiallelic positions ([52cc411](https://github.com/PharmGKB/PharmCAT/commit/52cc41157953110a4aa22349ee696bf31fc31de2)), closes [#102](https://github.com/PharmGKB/PharmCAT/issues/102)
* **preprocessor:** correct phasing status for multiallelic positions with mismatching alt ([8bbb9a6](https://github.com/PharmGKB/PharmCAT/commit/8bbb9a6f85bd1a446408e809c6f9dfaa54a79253)), closes [#102](https://github.com/PharmGKB/PharmCAT/issues/102)
* **preprocessor:** fix bug where VCF fields are not updated ([af604db](https://github.com/PharmGKB/PharmCAT/commit/af604db918a218c61a857f061d147eb206cd1a5e))
* **preprocessor:** harmonize arguments with PharmCAT ([b98c8ef](https://github.com/PharmGKB/PharmCAT/commit/b98c8ef05603526441cdb0ef65ad2371ab5b9464))
* **preprocessor:** keep the INFO/PX from refVcf in the output ([7bb3d00](https://github.com/PharmGKB/PharmCAT/commit/7bb3d0006e9aebd59ab61df48f4889500fd854f4))
* **preprocessor:** move bcftools and bgzip versions to global variables ([3b4ffda](https://github.com/PharmGKB/PharmCAT/commit/3b4ffda8de20e2b124fafd8a6556da22d71401f5))
* **preprocessor:** remove a line of test codes ([24a5881](https://github.com/PharmGKB/PharmCAT/commit/24a58810d9cb7db4c1780524621fd8fd5ddd6bbc))
* **preprocessor:** remove pre-existing .tbi index file ([7dd8c45](https://github.com/PharmGKB/PharmCAT/commit/7dd8c45ccf8f4893331592abce835870fbe470e1))
* **preprocessor:** revert the changes on .bgz which is not the default file extension of bgzip ([5b8e39c](https://github.com/PharmGKB/PharmCAT/commit/5b8e39c20c06b9cf70bccfcd46bbb388301abd06))
* **preprocessor:** sort non-PGx variants by genomic positions to prevent sorting error ([0b52234](https://github.com/PharmGKB/PharmCAT/commit/0b52234c4c4d4eee543a0ad92e93d7b2591589b7))
* **preprocessor:** update the default output suffix ([dc80699](https://github.com/PharmGKB/PharmCAT/commit/dc80699af3b37f1f8eeda4bb27283630d3b1ec94))
* **preprocessor:** use .bgz as a bgzipped file extension for clarify ([6b0e96f](https://github.com/PharmGKB/PharmCAT/commit/6b0e96f6ec6e4f4772b6a731416c5174c6231233))
* **preprocessor:** use sample IDs as the default output prefix ([b25b478](https://github.com/PharmGKB/PharmCAT/commit/b25b478e8f2ae5e3449e25b7a124b9eb255bc542))
* **preprocessor:** use the .bgz for bgzipped files ([9df0fa5](https://github.com/PharmGKB/PharmCAT/commit/9df0fa52681605da28898180469f79fff3672d6c))
* **preprocessor:** warn and quit if a lower version of bcftools is used ([49b3d29](https://github.com/PharmGKB/PharmCAT/commit/49b3d295f16d5e8e1e775056e64a23895e3202b5))
* **reporter:** add support for CPIC version ([b6b4714](https://github.com/PharmGKB/PharmCAT/commit/b6b4714f37615c389a545fde614851d02d732789))
* **reporter:** add validation for messages ([530c88b](https://github.com/PharmGKB/PharmCAT/commit/530c88b95b347101845f6daba343089f6fb3edfa))
* **reporter:** cleanup and use correct section name ([cea8f15](https://github.com/PharmGKB/PharmCAT/commit/cea8f15d5f9de66e5de532114c736a1fd35a444b))
* **reporter:** cleanup api ([87bf8c8](https://github.com/PharmGKB/PharmCAT/commit/87bf8c8a7b24c5bf3cbb7bf19743eaa59c03029e))
* **reporter:** correctly fill CYP2D6 copy number diplotype ([8d68b1e](https://github.com/PharmGKB/PharmCAT/commit/8d68b1ec49a3dfb089c77452118f8d6079f3fa91))
* **reporter:** don't store sample-specific info in DPWG data model ([06c3ed9](https://github.com/PharmGKB/PharmCAT/commit/06c3ed9f676e77194a5dbda1d6f4a2b9fae9dedf))
* **reporter:** extract TextConstants ([fdf51bb](https://github.com/PharmGKB/PharmCAT/commit/fdf51bb726a0659cd8a47a1d9303992742697c3c))
* **reporter:** fix allele presence matching for DPWG guidelines ([bdca135](https://github.com/PharmGKB/PharmCAT/commit/bdca135d95035dd9f840a1b82b34a38a0ce03a5e))
* **reporter:** fix bad references to ReportContext and DrugLink comparison ([c89ca70](https://github.com/PharmGKB/PharmCAT/commit/c89ca70d168bd48851bb0c993566a7afcaa2f4d9))
* **reporter:** fix bug with calculating possible genotypes for DPWG guidelines ([f0cc372](https://github.com/PharmGKB/PharmCAT/commit/f0cc372b7655236c1635e7cc906464b16f29477d))
* **reporter:** fix comparison of drug objects ([7289b0c](https://github.com/PharmGKB/PharmCAT/commit/7289b0c89fc5d2acdd9bbec7fc835c64c0337d79))
* **reporter:** fix duplication of drugs in gene summary ([d53fe0b](https://github.com/PharmGKB/PharmCAT/commit/d53fe0bb7a21337548c42aed3f465fdeea4e2450))
* **reporter:** fix function lookup when second allele is not present ([872f690](https://github.com/PharmGKB/PharmCAT/commit/872f690699167a4ac8685e09fca396e29a4ebb89))
* **reporter:** fix HLA group names to work with phenotype comparison ([3ab2356](https://github.com/PharmGKB/PharmCAT/commit/3ab23561ca44b9a2fdc2105d817de973bacd83b5))
* **reporter:** fix unclosed HTML element in template ([8b248ea](https://github.com/PharmGKB/PharmCAT/commit/8b248ea7b6cddce135983c936e29a79cff7eafb0))
* **reporter:** fix warfarin flowchart ([204bef9](https://github.com/PharmGKB/PharmCAT/commit/204bef92ceee40f88a1f65a9afd35d4e54433f35))
* **reporter:** force background colors in print media ([9d714bf](https://github.com/PharmGKB/PharmCAT/commit/9d714bf97c2eea217efec73992c67d5320af07c8))
* **reporter:** improve message handling to support source-specific messages and warfarin craziness ([9744883](https://github.com/PharmGKB/PharmCAT/commit/9744883f8c3c4eb6880d38b1d7b96165cf41cff8))
* **reporter:** improve outside call handling ([a580e32](https://github.com/PharmGKB/PharmCAT/commit/a580e322c6216de530a5ac07fe77523b3dadee47))
* **reporter:** make OutsideCallParser tolerate empty lines ([78c1b36](https://github.com/PharmGKB/PharmCAT/commit/78c1b36a4f246dd54407511e95dec650e04c1999))
* **reporter:** more fixes for problems introduced with Diplotype refactor ([300a775](https://github.com/PharmGKB/PharmCAT/commit/300a775fc10b919f33ca3316eaa572206def708b))
* **reporter:** more improvements and messaging tweaks ([07781a3](https://github.com/PharmGKB/PharmCAT/commit/07781a3a1763d7ff5d9d84e449e61b27a644696f))
* **reporter:** move DpydCaller into reporter.caller package ([8f7220e](https://github.com/PharmGKB/PharmCAT/commit/8f7220eb6ab7556bac2b32077e97f8d90ea88c78))
* **reporter:** move title, date, and version back into ReportContext ([ebb4aa3](https://github.com/PharmGKB/PharmCAT/commit/ebb4aa37d76fd176799a67f825cfd8efc72c7518))
* **reporter:** outside call now sets source diplotypes ([e2bc62e](https://github.com/PharmGKB/PharmCAT/commit/e2bc62e082a3518fac520e7aed6701dee1a1f46e))
* **reporter:** rename LeastFunctionUtils to DpydCaller ([1dabc8e](https://github.com/PharmGKB/PharmCAT/commit/1dabc8e0bfaae0c6863f51603372ff64c2db908d))
* **reporter:** render recommendation HTML properly ([1145479](https://github.com/PharmGKB/PharmCAT/commit/11454793f1f42fe6ce428eaa4b57d4c213e5494b))
* **reporter:** reorganize message handling ([06c4056](https://github.com/PharmGKB/PharmCAT/commit/06c4056582214029272d2db8c470e8ad42f27f82))
* **reporter:** reporter improvements ([f418167](https://github.com/PharmGKB/PharmCAT/commit/f418167028e79cd9c0697d20baa3427e20c8e5aa))
* **reporter:** show drug annotations in the report again ([e2daaf2](https://github.com/PharmGKB/PharmCAT/commit/e2daaf291841ed79664a3d581c450c11aeabfc60))
* **reporter:** standardize GSON usage and use Env instead of using new DefinitionReader ([31cf889](https://github.com/PharmGKB/PharmCAT/commit/31cf889eefc49f07e835a8de62b5406dde6d6c01))
* **reporter:** support inferring CYP2D6 copy number ([1448764](https://github.com/PharmGKB/PharmCAT/commit/1448764d62cc030d09206f88470d16bd02d778fd))
* **reporter:** update missing variant input display ([9e0fe5a](https://github.com/PharmGKB/PharmCAT/commit/9e0fe5a1700188b3e924bd65ceb105d78c4a3f22))
* **reporter:** use unicode GTE ([84843e3](https://github.com/PharmGKB/PharmCAT/commit/84843e32e0fe45230984e26d0009235b78ec057b))
* **reporter:** validate that PharmGKB annotations contain groups ([43340f8](https://github.com/PharmGKB/PharmCAT/commit/43340f8b27532c43401b5ae034f220b701a77684))
* **website:** fix typos and broken links ([5cdfc87](https://github.com/PharmGKB/PharmCAT/commit/5cdfc878157e574645bd41a4e9bac9f3709b691f))


### Performance Improvements

* **namedallelematcher:** cache allele lookups ([13731b2](https://github.com/PharmGKB/PharmCAT/commit/13731b2059a4f5a8228928e119c1bf682d2d0afd))
* **pharmcat:** skip reading/writing files if possible ([37c0491](https://github.com/PharmGKB/PharmCAT/commit/37c0491e3efbfccca9adbe67b095c673a6b01fe1))

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
