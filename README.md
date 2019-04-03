# CuratePlasmids

CuratePlasmids is a command-line tool for identifying bacterial plasmids from the [NCBI nucleotide](https://www.ncbi.nlm.nih.gov/nucleotide/) database, or from your own sequence assemblies.

# Table of contents

* [Introduction](#Introduction)
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Quick start](#Quick-start)
* [Background and methods](#Background-and-methods)
* [Options and usage](#Options-and-usage)
* [Output files](#Output-files)
* [FAQ](#faq)
* [Acknowledgements](#Acknowledgements)
* [License](#License)



# Introduction

Retrieving complete plasmid sequences from the NCBI nucleotide database requires quality-filtering to exclude partial plasmid sequences, or chromosomal sequences mis-annotated as plasmids. I previously outlined methods to curate NCBI plasmids ([Orlek _et al_. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28286183)). Following similar methods, CuratePlasmids allows users to first retrieve putative complete plasmid sequences from NCBI, and then characterise these sequences by detecting [replicon loci](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068535/) and [rMLST loci](https://pubmlst.org/rmlst/). As a result, genuine complete plasmid sequences can be identified. CuratePlasmids also allows users to characterise their in-house assembled sequences to distinguish plasmid and chromosomal contigs.<br>

[PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb) is an online database of curated plasmids, updated every ~3 months, so this is an easy way to get hold of NCBI plasmid sequences, and as an online database it comes with nice interactive features. However, there are reasons why you may want to use CuratePlasmids instead (see [FAQ](#faq) for details). Notably, CuratePlasmids is useful if you want to:
* Characterise in-house sequence data that has not yet been uploaded to NCBI.
* Retrive the most up-to-date set of plasmids from NCBI.
* Identify ['chromid'](https://www.ncbi.nlm.nih.gov/pubmed/20080407) sequences.
* Check the curation process at each step (which sequences are being excluded and why); and use this to fine-tune curation methods if desired. 


# Requirements


* Linux or MacOS (with the [Bash shell](https://en.wikibooks.org/wiki/Bash_Shell_Scripting#What_is_Bash?), which is the default shell on MacOS and many Linux distributions)
* [Python](https://www.python.org/) 3 is required for the `database_setup.py` executable (tested using Python 3.5); the `curate_plasmids.py` executable works with Python 2 (tested with Python 2.7) or Python 3
* [edirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
* The [rMLST database](https://pubmlst.org/rmlst/) (follow installation instructions below)

Note, the [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder_db) database for replicon typing is included in the repository (retrieved 23-Apr-2018). Both "enterobacteriaceae" and "gram_positive" components of the database are included. You can use a more recent version of PlasmidFinder if you wish (see [Options and usage](#Options-and-usage)).

# Installation

First install the repository:<br>

```bash
git clone https://github.com/AlexOrlek/ATCG.git
cd ATCG
```
You should find the executable scripts (`curateplasmids.py` and `database_setup.py`) within the repository directory. If you add the path of this directory to your [$PATH variable](https://www.computerhope.com/issues/ch001647.htm), then ATCG can be run by calling the executable scripts e.g. `curateplasmids.py [`*`arguments...`*`]` from any directory location.

__Installing the rMLST database__:

The ribsomal multi-locus sequence typing (rMLST) database comprises sequences of allelic ribosomal loci; these sequences can be used as chromosomal markers i.e to detect chromosomal as opposed to plasmid sequence. The database is free, but you will need to agree to an associated licence agreement which __forbids the distribution of the database__, in order to protect the Intellectual Property.

To gain access and install the database, follow these steps:
* Register for a PubMLST account if you do not already have one. The link to register is [here](https://pubmlst.org/bigsdb). Click on "Register for a site-wide account".
* Login to your account at https://pubmlst.org/bigsdb and request access to Ribosomal MLST genome and Ribosomal MLST locus/sequence definitions under Registrations. Additionally, email Keith Jolley (keith.jolley@zoo.ox.ac.uk) and request a consumer id and secret so that you'll be able to access the database programatically.
* Put your consumer id and secret into a text file (file name does not matter), with the id on the first line and the secret on the second. The file contents should look something like the below snippet:

```bash
efKXmqp2D0EBlMBkZaGC2lPf
F$M+fQ2AFFB2YBDfF9fpHF^qSWJdmmN%L4Fxf5Gur3
```

* Install the database by running the `database_setup.py` executable (with Python 3 set as your default Python version), providing the `-s` flag with the file containing the secret id and secret:<br>

```bash
database_setup.py -s secretfile.txt`
```

* Follow the instructions that appear on screeen. The rMLST database will be installed in the databases directory within the repository.



# Quick start

To retrive and curate bacterial plasmids from NCBI:

`curateplasmids.py -e first.last@email.com -o output-directory`


To retrieve and curate a custom set of NCBI accessions:

`curateplasmids.py -e first.last@email.com --accessions myaccessions.txt -o output-directory`


To curate your own plasmid sequences, provide an input multifasta file:

`curateplasmids.py --fasta samples.fasta -o output-directory`



# Options and usage

`curateplasmids.py --help` produces a summary of all the options.


`--enterobacdbpath` and `--gramposdbpath` flags can be provided with paths to your own PlasmidFinder enterobacteriaceae and gram_positive BLAST databases (created by running the `makeblastdb` command on the fasta files using [command line BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279688/)).<br>
The `--taxonomyquery` and `--datequery` flags allow the NCBI query to be customised according to the source organism of the accession and the date the accession was first added to NCBI. Note that the taxonomy id associated with a given organism can be found through the [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi). For example, the command below would retrieve _Klebsiella aerogenes_ plasmids, added since the start of 2017:<br>

`curateplasmids.py -e first.last@email.com -o output-directory -q '"Klebsiella aerogenes"[porgn:__txid548]' -d '"2017/01/01"[PDAT] : "3000"[PDAT]'`

The `-s` flag specifies which NCBI source database(s) to include; by default both refseq and genbank databases will be included (`refseq_genbank`) but refseq only can be specified (`refseq`).<br>
By default, the number of threads is 1, but multi-threading is recommended to reduce computing time (for BLAST searches); the number of threads to use is specified using the `-t` flag; the value must not exceed the number of threads available on your machine.<br>
The `--accessions` flag allows a user to bypass the NCBI query stage, and instead use a custom set of NCBI accessions.<br>
The `--retrieveaccessionsonly` flag outputs until `accessions_filtered.tsv` (see [Output file](output-files) and [Background and methods](#background-and-methods)).<br>
The `--retrievesequencesonly` flag outputs until `accessions_filtered_deduplicated.fa`, but does not run the more time-consuming BLAST-based filtering.<br>
As an example, if you wish to update an existing database with more recent accessions, you could run CuratePlasmids with the `--retrieveaccessionsonly` flag, and compare retrieved accessions with those in the existing database to identify novel putative plasmid accessions that you may wish to include. Then, you could run the next stage of CuratePlasmids by providing the set of novel putative plasmids to the `--accessions` flag to determine plasmid accessions to be included in the existing database.



# Output files

The below table shows the outputs from running the complete pipeline (including retrieval of accessions from NCBI) with default settings.

File/Directory               	       	  | Description
----------------------------------------- | --------------------------------------------------------------------------------------------- 
downloaddate.txt		    	  | a record of when the accessions were retrieved
accessions.tsv			    	  | putative plasmid accessions retrieved from NCBI
incompleteaccessions.tsv       	    	  | as above but not annotated as complete; these accessions are excluded
accessions_filtered.tsv        	    	  | accessions remaining, after filtering based on accession title text
excludedaccessions.tsv	     	    	  | accessions excluded, after filtering based on accession title text  
accessions_filtered_biosamples.tsv  	  | biosample accessions associated with accessions_filtered.tsv  
accessions_filtered_biosamplemetadata.tsv | biosample metadata (submitter name) associated with biosample accessions  
duplicateaccessions.tsv		    	  | accessions excluded, after removal of duplicate sequences with shared metadata 
accessions_filtered.fa              	  | sequences of accessions_filtered.tsv		    
accessions_filtered_deduplicated.fa 	  | sequences remaining after removal of duplicate sequences with shared metadata  
plasmidfinder/                	    	  | directory containing outputs from BLASTing remaining sequences against plasmid replicon loci
rmlst/                        	    	  | directory containing output from BLASTing remaining sequences against chromosomal rMLST loci
plasmids.fa		     	    	  | plasmid sequences  
plasmids.tsv		     	    	  | plasmid accessions
rmlstrepaccessions.tsv	     	    	  | accessions that have one or more replicon / rmlst loci detected (excluded from plasmids.fa)
rmlstonlyaccessions.tsv	     	    	  | accessions that have one or more rmlst loci (excluded from plasmids.fa)  

Accessions in rmlstrepaccessions.tsv are likely to be chromids. Accessions in rmlstonlyaccessions.tsv may be chromid sequences that don't contain a known plasmid replicon locus; alternatively, they may be chromosomal sequences mis-annotated as plasmids.


# Background and methods

For background information on curating plasmids see recent papers: [Orlek _et al._ (2017)](https://www.ncbi.nlm.nih.gov/pubmed/28286183) and [Galata _et al._ (2018)](https://academic.oup.com/nar/article/47/D1/D195/5149885). A brief outline of the steps of the complete pipeline (including retrieval from NCBI) is given below:

1. Putative complete plasmid accessions, along with accompanying information such as accession title are downloaded from NCBI nucleotide. To be considered a putative plasmid, the accession must be annotated as "plasmid" and "complete".
2. The accessions are filtered using a regular expression search of the accession title text. For example, titles including the words "gene", "transposon", or "synthetic vector" would be excluded.
3. The filtered accession sequences are downloaded as a FASTA file. If both Refseq and Genbank accessions have been retrieved, there will be duplicates; therefore, deduplication is conducted: accessions with identical sequences and shared metadata (same biosample accession id and/or submitter name) are deduplicated. When selecting a single accession from duplicates, Refseq accessions are favoured over Genbank accessions. 
4. The remaining sequences are BLASTed again the PlasmidFinder replicon database and the rMLST database.
5. If a sequence contains no rMLST loci then it is considered a plasmid and included in the plasmids.fa output file. Accessions with rMLST loci detected are recorded.

# FAQ

* **How does CuratePlasmids differ from previously published plasmid curation methods?**
I previously published a similar method for plasmid curation ([Orlek _et al._ 2017](https://www.ncbi.nlm.nih.gov/pubmed/28286183)), but compared with CuratePlasmids, the methods in the paper differ in several key ways:
    * I used MLST rather than rMLST to filter chromosomal accessions. MLST loci are more limited as a chromosomal marker since MLST schemes cover fewer taxa.
    * I used a less stringent and more convoluted approach to decide whether a plasmid sequence was compete: accessions were not required to be annotated as "complete" as long as the title text indicated a "complete" sequence; CuratePlasmids instead excludes any accession that is not annotated as "complete", but does not require explicit mention of completeness in the title text.<br>
    
    [Galata _et al._ (2018)](https://academic.oup.com/nar/article/47/D1/D195/5149885) created [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb) using methods similar to those of the CuratePlasmids pipeline. However, methods of Galata _et al._ differ somewhat; notably:
    * Although rMLST loci are used to filter chromosomal sequences, in contrast to methods of CuratePlasmids, sequences with up to 5 rMLST loci are included in the database. However, according to a recent review article ([di Cenzo & Finan 2017](https://mmbr.asm.org/content/81/3/e00019-17)), a plasmid-like sequence encoding one or more ribosomal loci should actually be considered a "chromid", which is biologically distinct from a plasmid ([Harrison _et al._ 2010](https://www.ncbi.nlm.nih.gov/pubmed/20080407)). CuratePlasmids makes a distinction between plasmids and chromids whereas PLSDB does not. 
* **Are there any caveats I should be aware of when using CuratePlasmids?**
CuratePlasmids relies on NCBI annotation for plasmid topology information (circular / linear). Ideally, a plasmid annotated as 'complete' and 'linear' should be a genuine linear plasmid, but the topology annotation should probably be treated cautiously; a 'linear' plasmid could represent a circular plasmid that failed to circularise after assembly.


# Acknowledgements

I am grateful to Dr Keith Jolley informing me about programmatic access to the rMLST database and for pointing me towards [ConFindr](https://olc-bioinformatics.github.io/ConFindr/) software, which implements programmatic access. The `database_setup.py` executable used in CuratePlasmids is based on a script from ConFindr.


# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)
