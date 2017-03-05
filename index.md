### Overview
[MendelBase](https://github.com/OpenMendel/MendelBase.jl) is a utility package that contains the base functions of [OpenMendel](https://openmendel.github.io). It provides the information needed to process keywords and read data from external files. MendelBase is necessary to run any OpenMendel analysis package, and you should first install MendelBase, before you install the analysis packages.

### Installation
*Note: The three OpenMendel packages (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [Search](https://openmendel.github.io/Search.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any other OpenMendel package will run. It is easiest if these three packages are installed in the above order and before any other OpenMendel package.*

Within Julia, use the package manager to install MendelBase:

    Pkg.clone("https://github.com/OpenMendel/MendelBase.jl.git")

This package supports Julia v0.4 and v0.5.

### Input Files
To make your genetic data available to OpenMendel, the various input files must be in a format
that OpenMendel can read. Listed here are all of OpenMendel's input file types. However, not every analysis option will use every type. (Example input files can be found in the docs subfolder of each of the analysis packages.)

* [Control File](#control-file): Specifies the names of your data input and output files and any optional parameters (*keywords*) for the analysis. (For a list of common keywords, see [Keywords Table](#keywords-table) below).
* [Pedigree File](#pedigree-file): Gives information about your individuals, such as name, parents' names, and sex. This may be a PLINK FAM format file.
* [Locus File](#locus-file): Names and describes the genetic loci in your data.
* [Phenotype File](#phenotype-file): Lists the relationships between the genotypes and the phenotypes at the loci in the locus file.
* [SNP Definition File](#snp-definition-file): Defines your SNPs with information such as SNP name, chromosome, position, allele names. This may be a PLINK BIM format file.
* [SNP Data File](#snp-data-file): Holds the genotypes for your data set. This must be a standard binary PLINK BED file in SNP-major format. If you have a SNP data file, you must also have a SNP definition file.

If an analysis option uses a SNP Data file (and thus also a SNP Definition file), then it does *not* use a Locus file (and thus it also does not use a Phenotype file). Similarly, if an analysis option uses a Locus file, then it does not use a SNP Data nor SNP Definition file.

### Control files<a id="control-file"></a>
The Control file is a text file consisting of keywords and their assigned values. The format for each line of the Control file is:

	Keyword = Keyword_Value(s)

The names of keywords are *not* case sensitive, but the keyword values *may* be case sensitive. Below is an example Control file (any text after "#" is treated as comments and ignored):

	# Input and Output files.
	#
	locus_file = gamete competition LocusFrame.txt
	pedigree_file = gamete competition PedigreeFrame.txt
	output_file = gamete competition Output.txt
	#
	# Analysis parameters for the Gamete Competition option.
	#
	trait = ACE
	affected_designator = 1
	standard_errors = true

In the example above, there are six keywords. The first three specify input and output filenames: *locus_file*, *pedigree_file*, and *output_file*. The last three specify analysis parameters: *trait*,  *affected_designator*, *standard_errors*. The text after the "=" are the keyword values, for example, *ACE*, *1*, and *true*. For keywords that can take a set of values, the values are separated by commas.

### Keywords<a id="keywords-table"></a>
This is a list of the general OpenMendel keywords common to most analysis package. The keywords relevant to a specific analysis option are covered in detail in the documentation for that package. Again, the names of keywords are *not* case sensitive, but the keyword values *may* be case sensitive. (Double quotes (") denote strings; single quotes (') denote a single character.)

 Keyword          |   Default Value    | Allowed Values |  Short Description       
----------------  |  ----------------  |  ----------------  |  ----------------
affected_designator  |  affected  |    |  Value used to designate affecteds at the trait field
allele_pseudo_count  |  0.1  |  Real  |  Count of pseudo alleles or haplotypes to use as a prior in a Bayesian method for estimating frequencies
allele_separator  |  /  |  1 Character  |  Separates alleles in unordered genotypes
complexity_threshold  |  5.00E+07  |  Real  |  Threshold at which pedigrees with greater complexity are not analyzed because the computation would take too long
eliminate_genotypes  |  FALSE  |  TRUE, FALSE  |  Indicator of whether genotypes should be eliminated
female  |  "female", "f", 'f', "2", '2', 2  |   Female designators separated by commas  |  Female designators
field_separator  |  ','  |  1 Character  |  Character used to separate values in data files
genetic_map_function  |  Haldane   |   Haldane or Kosambi    |  Function for converting Morgans to recombination fractions
locus_file  |    |  Name of existing Locus file  |  Name of input Locus file
lump_alleles  |  FALSE  |  TRUE, FALSE  |  Indicator of whether alleles should be combined
male  |  "male", "m", 'm', "1", '1', 1  |  Male designators separated by commas  |  Male designators
new_pedigree_file  |    |  Filename  |  Name of file that will hold OpenMendel generated pedigree data
ordered_allele_separator  |  \| (vertical bar) |  1 Character  |  Separates alleles in ordered genotypes
output_file  |  Mendel_Output.txt  |  Filename  |  Name of file that will hold OpenMendel analysis results
pedigree_file  |    |  Name of existing Pedigree file  |  Name of input Pedigree file
phenotype_file  |    |  Name of existing Phenotype file  |  Name of input Phenotype file
plink_input_basename  |    |  Basename of existing .bed, .bim, and .fam files  |  Basename of PLINK input files
plink_output_basename |    |    |  Basename of PLINK output files
populations  |    |  Names of populations separated by commas  |  Names of populations
product_mode  |  TRUE  |  TRUE, FALSE  |  
seed  |  1234  |  Integer  |  Seed value for random number generator
snpdata_file  |    |  Name of existing SNP Data file (in .bed format)  |  Name of input SNP Data file
snpdefinition_file  |    |  Name of existing SNP Definition file  |  Name of input SNP Definition file
standard_errors  |  FALSE  | TRUE, FALSE   |  Indicator whether standard errors should be output
trait  |    |    |  Name of trait column in pedigree file (use Trait for .fam files)

### General Data File Format
All data files are text files with one exception: if used, the SNP data file should be a binary file in PLINK BED file format. OpenMendel will accept [PLINK format](http://zzz.bwh.harvard.edu/plink) FAM and BIM files. All non-PLINK OpenMendel data files list their data in columns, with a header row with names for each column. The columns need not be in any fixed order. Many column have required, case sensitive names that need to be in the header row. The column separators can be any single character, specified in the control file with the keyword *field_separator*. The default separator is a *comma*. If a data element contains the field separator (e.g., a comma) it must be enclosed in double quotes, or it will be read as multiple fields (e.g., "A/A, A/O"). Missing data is represented by a blank or the two characters NA. In PLINK files the missing data symbols are 0 or -9.
	
### Pedigree File<a id="pedigree-file"></a>
A Pedigree file is always required. OpenMendel will accept [PLINK format](http://zzz.bwh.harvard.edu/plink) FAM pedigree files, in which case the file must have the extension *.fam*. The OpenMendel pedigree file format has a header row with some of the following field names (in any order): Pedigree, Individual (or Person), Mother, Father, Sex, and Twin. The only required column is the one labeled either *Individual* or *Person* that provides the samples names.

If the Pedigree column is missing, all rows are treated as separate pedigrees of size one, and all individuals are unrelated. Within each pedigree, every individual must have a unique name. If the Mother and Father columns are missing, all individuals are treated as founders.  If the Mother and Father columns exist, each individual must *either* have both a named mother and father, *or* both parental name fields must be blank. If both fields are blank, that individual is treated as a founder. If the Sex field is missing, all individuals are treated as female. If a Sex column exists but the sex is missing for an individual, that individual is treated as female. If the Twin field is missing there are no monozygotic (MZ) twins in the data set. If the Twin field is present and a set of individuals in the same pedigree all have the same non-blank value in their Twin field, they are all considered MZ siblings.

If the analysis uses a [Locus file](#locus-file) to describe the markers, the marker genotypes are contained in the Pedigree file, and the marker names in the Locus and Pedigree files must agree. If the analysis uses a [SNP Definition file](#snp-definition-file) to describe the markers, the SNP genotypes are contained in the [SNP Data file](#snp-data-file), not the Pedigree file. The Pedigree file may also include optional columns such as traits and [population](#populations) ancestries.

**Example Pedigree File**

	Pedigree,Person,Mother,Father,Sex,Twin,ABO,HtMeters
	Clinton,Bill,,,male,0,A,1.88
	Clinton,Hillary,,,female,0,O,NA
	Clinton,Chelsea,Hillary,Bill,female,0,A,NA
	Obama,Barack,,,male,0,O,1.85
	Obama,Michelle,,,female,0,B,NA
	Obama,Malia,Michelle,Barack,female,0,O,NA
	Obama,Sasha,Michelle,Barack,female,0,B,NA

### Locus File<a id="locus-file">
The Locus file defines your genetic markers if you are working with a relatively small number of markers. (If you have a dense, genome-wide set of SNPs, then you should use a [SNP Definition file](#snp-definition-file) to define your SNPs.) The Locus file must have a header row with the field names (in any order): *Locus*, *Allele*, *Chromosome*, and *Frequency*. You may also include locus position information with the optional fields *Morgans*, *FemaleMorgans* and *MaleMorgans*, and *Basepairs*. The Locus file contains one line for each allele at each locus. The lines for all the alleles for a particular locus should be contiguous.

**Example Locus File**

	Locus,Allele,Chromosome,Frequency
	ABO,A,9,0.64,0.27
	ABO,B,9,0.34,0.06
	ABO,O,9,1.35,NA
	Rh,D,1,2.42,0.61
	Rh,d,1,0.58,0.39
	SNP,1,9,0,NA
	SNP,2,9,0,NA

### Phenotype File<a id="phenotype-file">

If there are loci in the Locus file that are non-codominant, then the relationships between the genotypes and the phenotypes for these loci are listed in the Phenotype File. The Phenotype File is in the standard table format with a header row. The header row must have the required field names *Locus*, *Phenotype*, and *Genotypes* in any order. Each line names a locus and a phenotype used at that locus. Also listed on that line are all the genotypes consistent with that phenotype, in a quoted list. For example, the line `ABO,A,"A/A,A/O"` indicates the genotypes A/A and A/O are both valid for the phenotype A at the ABO locus.

**Example Phenotype File**

	Locus,Phenotype,Genotypes
	ABO,A,"A/A,A/O"
	ABO,B,"B/B,B/O"
	ABO,AB,"A/B"
	ABO,O,"O/O"
	Rh,+,"D/D,D/d"
	Rh,-,"d/d"
	Xg,+,"+/+,+/-"
	Xg,-,"-/-"

### SNP Definition File<a id="snp-definition-file">
The SNP Definition file defines the SNPs with information such as SNP name, chromosome, position, and allele names. OpenMendel will accept [PLINK format](http://zzz.bwh.harvard.edu/plink) BIM files, in which case the file must have the extension *.bim*. The OpenMendel SNP Definition file format has a header row (in any order) with some of the field names: *SNP* or *Locus*, *Chromosome*, *CentiMorgans*, *Basepairs*, *Allele1*, and *Allele2* (these last two list the allele names). If neither a SNP nor Locus field is present, then the n SNPs are named 1, 2, 3, ..., n. If the Chromosome field is missing altogether or contains a missing value at some SNP, it is assumed to be autosomal. If present, the Chromosome field should have a value from 1 to 22 or the value X. If the CentiMorgan or Basepair fields are missing, they take the value 0. If Allele1 is missing it is assigned the value "1"; a missing Allele2 becomes "2".

**Example SNP Definition File**

	Chromosome,SNP,CentiMorgans,Basepairs,Allele1,Allele2
	X,rs311165,0,2699555,C,A
	X,rs28579419,0,2699645,T,G
	X,rs2306736,0,2700027,C,T
	X,rs5939320,0,2700202,A,G
	X,rs4892819,0,2700613,G,A
	X,rs5982853,0,2702143,C,T
	X,rs73433431,0,2702698,T,C

### SNP Data File<a id="snp-data-file">
The SNP Data file holds the SNP genotypes. The SNP Data file must be a standard [PLINK format](http://zzz.bwh.harvard.edu/plink) binary BED file in SNP-major format. If you have a SNP data file, you must have a SNP definition file, but it need not be in PLINK format. (We do not list here an example SNP Data file because they are binary files, and thus not human readable.)

### PLINK Basename
You may use the keyword *plink_input_basename* in the Control file to specify a name for your PLINK format files. If you use this keyword you do not need to name each PLINK file in the Control file, but you must supply each of the three relevant PLINK files, using that basename:  *your_basename.bed* (SNP data file), *your_basename.bim* (SNP definition file), and *your_basename.fam* (pedigree file).

### Populations<a id="populations">
If your individuals are of known mixed ethnicity, OpenMendel will allow you to adjust for this in your analysis by using the *populations* keyword. The relevant populations should be named in your [Control file](#control-file), fractional ancestries should be assigned to the individuals in the [Pedigree file](#pedigree-file), and population-specific allele frequencies should be assigned to the markers in the [Locus file](#locus-file).

For example, the [Control file](#control-file) might be:

	# Input and Output files.
	#
	locus_file = genedropping LocusFrame.txt
	pedigree_file = genedropping PedigreeFrame.txt
	phenotype_file = genedropping PhenotypeFrame.txt
	new_pedigree_file = genedropping NewPedigreeFrame.txt
	#
	# Analysis parameters for the Gene Dropping option.
	#
	repetitions = 2
	populations = European, African, Chinese
	
The [Pedigree file](#pedigree-file) should have a column for each population, with fractional ancestries assigned at each individual. The fractional ancestries need not sum to 1, but will be normalized to 1 in the analysis. An example [Pedigree file](#pedigree-file) is:

	Pedigree,Person,Mother,Father,Sex,Twin,European,African,Chinese,ABO
	Clinton,Bill,,,male,0,1,0,0,A
	Clinton,Hillary,,,female,0,1,0,0,O
	Clinton,Chelsea,Hillary,Bill,female,0,NA,NA,NA,A
	Obama,Barack,,,male,0,0.5,0.5,0,O
	Obama,Michelle,,,female,0,0.3,0.7,0,B
	Obama,Malia,Michelle,Barack,female,0,NA,NA,NA,O
	Obama,Sasha,Michelle,Barack,female,0,NA,NA,NA,B

The [Locus file](#locus-file) should have, instead of one column labeled *Frequency*, a frequency column for each population, where the header for the column is the name of that population as listed in the Control file. Each row then lists for that allele its frequency in each population. The population frequencies need not sum to 1, but will be normalized to 1 in the analysis. An example [Locus file](#locus-file) is:

	Locus,Allele,Chromosome,European,African,Chinese
	ABO,A,9,0.27,0.18,0.19
	ABO,B,9,0.06,0.11,0.17
	ABO,O,9,NA,0.71,0.64

### Output Files
All OpenMendel analysis options will create a results file. If no name for this file is specified in the [Control file](#control-file), it will be given the default name *Mendel_Output.txt*. A number of OpenMendel analysis options can create additional output files, such as a plot file or a new pedigree data file, based on the analysis. Newly created data files can be used as input files in subsequent analyses. To create these additional output files, the filename is specified by a keyword in the Control file. For example:

	# Input and Output files.
	#
	plink_input_basename = gwas 1 data
	output_file = gwas 1 Output.txt
	manhattan_plot_file = ManPlotOutput_1.png
	new_pedigree_file = NewPedigreeFrame.txt
	#
	# Analysis parameters for the GWAS option.
	#
	regression = linear
	regression_formula = Trait ~ Sex

### Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ### Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

### Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
