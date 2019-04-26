# bacterialBercow

bacterialBercow is a command-line tool for retrieving complete bacterial plasmids from the [NCBI nucleotide](https://www.ncbi.nlm.nih.gov/nucleotide/) database; and for characterising contigs from your own bacterial sequence assemblies.

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

__Retrieving complete bacterial plasmids from NCBI__:<br>
Retrieving complete plasmid sequences from the NCBI nucleotide database requires quality-filtering to exclude partial plasmid sequences, or chromosomal sequences mis-annotated as plasmids. I previously outlined methods to curate NCBI bacterial plasmids ([Orlek _et al_. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28286183)). Following similar methods, bacterialBercow allows users to first retrieve putative complete bacterial plasmid sequences from NCBI, and then characterise these sequences by detecting plasmid [replicon loci](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068535/) and ribosomal multi-locus sequence typing ([rMLST](https://pubmlst.org/rmlst/)) loci. As a result, genuine complete plasmid sequences can be identified.

__Characterising in-house bacterial sequence assemblies__:<br>
bacterialBercow also allows users to characterise their own in-house assembled sequences, using plasmid replicon typing and rMLST typing. This can help in distinguising plasmid and chromosomal sequences; in addition, rMLST can be used to determine bacterial species.<br>

__Why use bacterialBercow__?<br>
[PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb) is an online database of curated plasmids, updated every ~3 months, so this is an easy way to get hold of NCBI plasmid sequences, and as an online database it comes with nice interactive features. However, there are reasons why you may want to use bacterialBercow instead (see [FAQ](#faq) for details). Notably, bacterialBercow is useful if you want to:
* Characterise in-house sequences that have not yet been uploaded to NCBI, using replicon typing and rMLST typing. Note that in addition to rMLST loci being useful markers of non-plasmid (chromosomal or [chromid](https://www.ncbi.nlm.nih.gov/pubmed/20080407)) sequence, the rMLST typing scheme can be used to determine bacterial species, given a complete chromosome (and it may be possible to determine species or genus-level taxonomic information from incomplete chromosomal sequence, if sufficient rMLST loci are represented). As far as I know, there are currently no other command-line tools for rMLST-based taxonomy (only the [rMLST website](https://pubmlst.org/rmlst/)). 
* Retrieve the most up-to-date set of plasmids from NCBI.
* Check the curation process at each step (which sequences are being excluded and why); and use this to fine-tune curation methods if desired. 


# Requirements


* Linux or MacOS (with the [Bash shell](https://en.wikibooks.org/wiki/Bash_Shell_Scripting#What_is_Bash?), which is the default shell on MacOS and many Linux distributions)
* [Python](https://www.python.org/) 3 is required for the `database_setup.py` executable (tested using Python 3.5); the `order.py` executable works with Python 2 (tested with Python 2.7) or Python 3
* [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK52640/#_chapter1_Installation_) (tested using version 2.6.0)
* [edirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/)
* A computer with internet access via the HTTPS protocol - required for retrieving data from NCBI.
* [bioawk](https://github.com/lh3/bioawk)
* The [rMLST database](https://pubmlst.org/rmlst/) (follow installation instructions below)

Note, the [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder_db) database for replicon typing is included in the repository (retrieved 9th-Apr-2019). Both "enterobacteriaceae" and "gram_positive" components of the database are included. You can use a more recent version of PlasmidFinder if you wish (see [Options and usage](#Options-and-usage)).

# Installation

First install the repository:<br>

```bash
git clone https://github.com/AlexOrlek/bacterialBercow.git
cd bacterialBercow
```
You should find the executable scripts (`order.py` and `database_setup.py`) within the repository directory. If you add the path of this directory to your [$PATH variable](https://www.computerhope.com/issues/ch001647.htm), then bacterialBercow can be run by calling the executable scripts e.g. `order.py [`*`arguments...`*`]` from any directory location.

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
database_setup.py -s secretfile.txt
```

* Follow the instructions that appear on screeen. The rMLST database will be installed in the databases directory within the repository.



# Quick start

To retrieve and curate bacterial plasmids from NCBI:

`order.py -e first.last@email.com -o output-directory`


To retrieve and curate a custom set of NCBI accessions:

`order.py -e first.last@email.com --accessions accessions.txt -o output-directory`


To characterise your own bacterial sequences, provide an input multi-FASTA file:

`order.py --inhousesequences samples.fasta -o output-directory`



# Options and usage

`order.py --help` produces a summary of all the options. All commands must include the `-o` flag specifying output directory. The other key flags are summarised below:

__NCBI query and retrieval options__:<br>
A contact email address must be provided using the `-e` flag; this is so that NCBI have someone to notify in case the software causes problems for their server.<br>
The `--taxonomyquery` and `--datequery` flags allow the NCBI query to be customised according to the source organism of the accession and the date the accession was first added to NCBI. Note that the taxonomy id associated with a given organism can be found through the [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi). For example, the command below would retrieve _Klebsiella aerogenes_ plasmids, added since the start of 2017:<br>

`order.py -e first.last@email.com -o output-directory -taxonomyquery '"Klebsiella aerogenes"[porgn:__txid548]' -datequery '"2017/01/01"[PDAT] : "3000"[PDAT]'`

The `-s` flag specifies which NCBI source database(s) to include; by default both Refseq and Genbank databases will be included (`refseq_genbank`) but Refseq only can be specified (`refseq`).<br>
The `--deduplicationmethod` flag specifies how identical sequences should be deduplicated (see [Background and methods](#background-and-methods) for details).<br>


__Replicon typing and rMLST typing options__:<br>
By default, the number of threads is 1, but multi-threading is recommended to reduce the computing time of the BLAST searches; the number of threads to use is specified using the `-t` flag; the value must not exceed the number of threads available on your machine.<br>
The `--typing` flag is applicable if the `--inhousesequences` flag is provided. By default both replicon and rMLST typing will be conducted on in-house sequences, but a user can specify only `replicon` typing or only `rmlst` typing. For accessions retrieved from NCBI, both replicon and rMLST typing are performed since this helps to distinguish plasmids from chromosomal sequence and chromids.<br>
`--enterobacdbpath` and `--gramposdbpath` flags can be provided with paths to your own PlasmidFinder enterobacteriaceae and gram_positive BLAST databases (created by running the `makeblastdb` command on the FASTA files using [command line BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279688/)).<br>


__Pipeline step options specifying starting and stopping points__:<br>
If provided, the `--retrieveaccessionsonly` flag will stop the pipeline after `accessions_filtered.tsv` is produced (see [Output file](output-files) and [Background and methods](#background-and-methods)).<br>
If provided, the `--retrievesequencesonly` flag will stop the pipeline after `accessions_filtered_deduplicated.fa` is produced.<br>
If provided, the `--restartwithsequences` flag will re-start the pipeline from the point where `--retrievesequencesonly` stopped the pipeline. The output directory specified must be the same as the output directory that was previously specified when the pipeline was run with the `--retrievesequencesonly` flag.<br>
The `--accessions` flag allows a user to bypass the NCBI query stage, and instead use a custom set of NCBI accessions.<br>
The `--inhousesequences` flag allows a user to provide their own multi-FASTA file of sequences which will be characterised using replicon typing and rMLST typing.<br>
<br>
_Example 1: you wish to update an existing database with more recent accessions:_<br>
You could run bacterialBercow with the `--retrieveaccessionsonly` flag, and compare retrieved accessions with those in the existing database to identify novel putative plasmid accessions that you may wish to include. Then, you could run the next stage of bacterialBercow by providing the set of novel putative plasmids to the `--accessions` flag to determine plasmid accessions to be included in the existing database.<br>
_Example 2: you have access to a computer cluster which can be run with lots of `--threads` but for security reasons, the HTTPS protocol (required for accessing data from NCBI) is not permitted:_<br>
You could run bacterialBercow in two stages. First, on your own computer, run bacterialBercow with the `--retrievesequencesonly` flag. Then, with the same output directory (`-o` flag) specified, run the typing stage of the pipeline on your computer cluster by providing the `--restartwithsequences` flag along with lots of `--threads`.



# Output files

The below table shows the outputs from running the complete pipeline (including retrieval of accessions from NCBI) with default settings.

File/Directory            	    | Description
----------------------------------- | --------------------------------------------------------------------------------------------- 
downloaddate.txt		    | a record of when the accessions were retrieved
accessions.tsv			    | putative plasmid accessions retrieved from NCBI
incompleteaccessions.tsv       	    | as above but not annotated as complete; these accessions are excluded
accessions_filtered.tsv        	    | accessions remaining, after filtering based on accession title text
excludedaccessions.tsv	     	    | accessions excluded, after filtering based on accession title text  
accessions_filtered_dblinks.tsv     | BioSample and BioProject accession ids associated with nucleotide accessions from accessions_filtered.tsv  
accessions_filtered_metadata.tsv    | BioSample metadata (submitter name and owner name) associated with BioSample accessions
accessions_filtered_refseq_gb.tsv   | Refseq accessions and their cognate Genbank accessions
identicalaccessions.tsv		    | accessions with identical sequences; corresponding submitter metadata and BioProject accession ids 
accessions_filtered.fa              | sequences of accessions_filtered.tsv		    
accessions_filtered_deduplicated.fa | sequences remaining after removal of duplicate sequences  
plasmidfinder/                	    | directory containing outputs from BLASTing remaining sequences against plasmid replicon loci
rmlst/                        	    | directory containing output from BLASTing remaining sequences against chromosomal rMLST loci
plasmids.fa		     	    | plasmid sequences  
plasmids.tsv		     	    | information on plasmid sequences including replicon typing
nonplasmids.tsv			    | information on non-plasmid sequences including replicon typing and rMLST typing

The below table shows outputs from running the pipeline with the `--inhousesequences` flag provided.

File/Directory            	    | Description
----------------------------------- | --------------------------------------------------------------------------------------------- 
plasmidfinder/                	    | directory containing outputs from BLASTing sequences against plasmid replicon loci
rmlst/                        	    | directory containing output from BLASTing sequences against chromosomal rMLST loci
seqlengths.tsv			    | sequence names (from FASTA headers) and corresponding sequence lengths
typing.tsv			    | information on sequences including replicon typing and rMLST typing



# Background and methods

For background information on curating NCBI plasmids see recent papers: [Orlek _et al._ (2017)](https://www.ncbi.nlm.nih.gov/pubmed/28286183) and [Galata _et al._ (2018)](https://academic.oup.com/nar/article/47/D1/D195/5149885). A brief outline of the steps of the complete pipeline (including retrieval from NCBI) is given below:

1. Putative complete plasmid accessions, along with accompanying information such as accession title are downloaded from NCBI nucleotide. To be considered a putative plasmid, the accession must be annotated as "plasmid" and "complete".
2. To select for complete plasmid genomes, the accessions are filtered using a regular expression search of the accession title text. For example, titles including the words "gene", "transposon", or "synthetic vector" would be excluded.
3. The filtered accession sequences are downloaded as a FASTA file. If both Refseq and Genbank accessions have been retrieved, there will be duplicates; therefore, deduplication is conducted, leaving only one sequence from the duplicate set of sequences with identical nucleotide sequence, favouring retention of Refseq over Genbank accessions. More specifically, deduplication is conducted according to the argument provided to the `--deduplicationmethod` flag:
    * `all`: all duplicate sequences are filtered, leaving only one representative sequence.
    * `bioproject`: sequences are deduplicated if they share the same BioSample accession id, and the same BioProject accession id (default method).
    * `submitter`: sequences are deduplicated if they share the same BioSample accession id, and the same submitter metadata (submitter contact name and submitter affiliation name, from the [BioSample "Owner"](https://www.ncbi.nlm.nih.gov/books/NBK169436/) field)
    
    The aim of the `bioproject` and `submitter` methods is to allow for interesting duplicates to be retained - identical plasmids from distinct samples and sequencing projects, that may represent transmission of a conserved plasmid. Nucleotide accessions originating from the same biological sample should share the same BioSample accession; nucleotide accessions originating from the same sequencing project should share the same primary BioProject id ([Pruitt 2011](https://www.ncbi.nlm.nih.gov/books/NBK54015/), [Barrett 2012](https://www.ncbi.nlm.nih.gov/pubmed/22139929), [Benson 2013](https://academic.oup.com/nar/article/41/D1/D36/1068219)). However, in Refseq, higher-level primary BioProject accessions [can be created by NCBI staff](https://www.ncbi.nlm.nih.gov/bioproject/docs/faq/#what-is-project-type), e.g. [PRJNA224116](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA224116). This breaks the correspondence between Refseq BioProject accession and sequencing project; therefore, for Refseq accessions the original BioProject accession id is retrieved via the cognate Genbank accession.    
4. The remaining sequences are BLASTed again the PlasmidFinder replicon database and the rMLST database.
5. If a sequence contains no rMLST loci then it is considered a plasmid and included in the plasmids.fa output file (although see the [FAQ](faq) for the potential limitations of this approach). Accessions with rMLST loci detected are recorded; these could be chromid or chromosomal sequences.

# FAQ

* **Why is the tool called bacterialBercow?**
At the time of writing, [John Bercow](https://en.wikipedia.org/wiki/John_Bercow) is the [Speaker of the House of Commons](https://en.wikipedia.org/wiki/Speaker_of_the_House_of_Commons_(United_Kingdom)); he is responsible for bringing order to UK parliamentary debates in [unruly times](https://www.youtube.com/watch?v=EY7EIZl4raY), and does so by yelling "Order! Ordeerr!". As I explained in my "Ordering the mob" paper, these are unruly times for microbiologists too - faced with a deluge of bacterial genome sequences, available via NCBI and from in-house sequencing projects. 
* **How do the methods of bacterialBercow differ from previously published methods for retrieving and curating NCBI plasmids?**
I previously published a similar method for plasmid curation ([Orlek _et al._ 2017](https://www.ncbi.nlm.nih.gov/pubmed/28286183)), but compared with bacterialBercow, the methods in the paper differ in several key ways:
    * I used [MLST](https://pubmlst.org/general.shtml) rather than rMLST to filter chromosomal accessions. MLST loci are more limited as a chromosomal marker since MLST schemes cover fewer taxa.
    * I used a less stringent and more convoluted approach to decide whether a plasmid sequence was compete: accessions were not required to be annotated as "complete" as long as the title text indicated a "complete" sequence; bacterialBercow instead excludes any accession that is not annotated as "complete", but does not require explicit mention of completeness in the title text.<br>
    
    [Galata _et al._ (2018)](https://academic.oup.com/nar/article/47/D1/D195/5149885) created [PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb) using methods similar to those of the bacterialBercow pipeline. However, methods of Galata _et al._ differ somewhat; notably:
    * Although rMLST loci are used to filter chromosomal sequences, in contrast to methods of bacterialBercow, sequences with up to 5 rMLST loci are included in the database. However, according to a recent review article ([di Cenzo & Finan 2017](https://mmbr.asm.org/content/81/3/e00019-17)), a plasmid-like sequence encoding one or more ribosomal loci should actually be considered a "chromid", which is biologically distinct from a plasmid ([Harrison _et al._ 2010](https://www.ncbi.nlm.nih.gov/pubmed/20080407)). bacterialBercow makes a distinction between plasmids and chromids whereas PLSDB does not.
    * PLSDB excludes all duplicate sequences (any sequences with a [mash distance](https://mash.readthedocs.io/en/latest/index.html) of 0). By default, bacterialBercow excludes identical sequences except those derived from a different BioSample and BioProject. The latter are likely to represent interesting cases of transmission of short conserved plasmids. Information on all duplicates is recorded.
* **Are there any caveats I should be aware of when using bacterialBercow?**
    * To curate NCBI plasmids, bacterialBercow relies on NCBI annotations for plasmid topology information (circular / linear). Ideally, a plasmid annotated as 'complete' and 'linear' should be a genuine linear plasmid, but the topology annotation should probably be treated cautiously, particularly since linear is the [default value when submitting to NCBI](https://www.ncbi.nlm.nih.gov/books/NBK293904/); also, a 'linear' plasmid could represent a circular plasmid that failed to circularise after assembly.
    * Characterising sequences using plasmid replicon typing and rMLST can help towards their categorisation (as plasmid / chromosome / chromid). For example, a sequence encoding >50 rMLST loci is likely to be a complete or near-complete chromosome ([Jolley _et al._ 2012](https://mic.microbiologyresearch.org/content/journal/micro/10.1099/mic.0.055459-0#tab2)). Likewise, any sequence with one or more rMLST loci is either a chromid or chromosomal sequence. However, the typing results should be interpreted cautiously. A chromid is a plasmid-like sequence that encodes some essential bacterial genes ([di Cenzo & Finan 2017](https://mmbr.asm.org/content/81/3/e00019-17)). rMLST loci (which encode ribosomal proteins) are not the only group of essential genes ([Gil _et al._ 2004](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC515251/)); therefore, absence of rMLST loci does not necessarily mean that a sequence is a plasmid - it could be a chromid sequence encoding other essential genes. Future versions of bacterialBercow could incorporate more comprehensive sets of bacterial core genes in order to better distinguish chromid sequences.


# Acknowledgements

I am grateful to Dr Keith Jolley for informing me about programmatic access to the rMLST database ([Jolley _et al._ 2017](https://academic.oup.com/database/article/doi/10.1093/database/bax060/4079979)) and for pointing me towards [ConFindr](https://olc-bioinformatics.github.io/ConFindr/) software, which implements programmatic access. The `database_setup.py` executable used in bacterialBercow is based on a script from ConFindr.


# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)
