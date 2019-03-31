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

Retrieving plasmid sequences from the NCBI nucleotide database requires quality-filtering to exclude partial plasmid sequences, or chromosomal sequences mis-annotated as plasmids. I previously outlined methods to curate NCBI plasmids ([Orlek _et al_. 2017](https://www.ncbi.nlm.nih.gov/pubmed/28286183)). Following similar methods, CuratePlasmids allows users to retrieve putative plasmid sequences from NCBI and then characterise them by detecting [replicon loci](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068535/) and [rMLST loci](https://pubmlst.org/rmlst/), so that plasmid sequences can be distinguished. CuratePlasmids also allows users to characterise their in-house assembled sequences.<br>

[PLSDB](https://ccb-microbe.cs.uni-saarland.de/plsdb) is an online database of curated plasmids, updated every ~3 months, so this is an easy way to get hold of NCBI plasmid sequences. However, there are reasons why you may want to use CuratePlasmids instead (see [FAQ](#faq) for details). Notably, CuratePlasmids is useful if you want to:
* Characterise in-house sequence data that has not yet been uploaded to NCBI.
* Retrive the most up-to-date set of plasmids from NCBI.
* Identify ['chromid'](https://www.ncbi.nlm.nih.gov/pubmed/20080407) sequences.
* Check the curation process at each step (which sequences are being excluded and why); and use this to fine-tune curation methods if desired. 


# Requirements


* Linux or MacOS (with the [Bash shell](https://en.wikibooks.org/wiki/Bash_Shell_Scripting#What_is_Bash?), which is the default shell on MacOS and many Linux distributions)
* [Python](https://www.python.org/) 3 (tested using Python 3.5)
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


`--enterobacdbpath` and `--gramposdbpath` flags can be provided with paths to your own PlasmidFinder enterobacteriaceae and gram_positive BLAST databases (created by running the `makeblastdb` command on the fasta files using [command line BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279688/))<br>
The `--taxonomyquery` and `--datequery` flags allow the NCBI query to be customised according to the source organism of the accession and the date the accession was first added to NCBI. Note that the taxonomy id associated with given organism can be found through the [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi). For example, the command below would retrieve _Klebsiella aerogenes_ plasmids, added since the start of 2017:<br>

`curateplasmids.py -e first.last@email.com -o output-directory -q '"Klebsiella aerogenes"[porgn:__txid548]' -d '"2017/01/01"[PDAT] : "3000"[PDAT]'`

The `--accessions` flag allows a user to bypass the NCBI query stage, and instead use a custom set of NCBI accessions.<br>
The `--retrieveaccessionsonly` flag retrives accessions and runs the initial title text-based filtering, but not the more time-consuming BLAST-based filtering (see [Background and methods](#background-and-methods)).<br>
Therefore, if you wish to update an existing database with more recent accessions, you could run CuratePlasmids with the `--retrieveaccessionsonly` flag, and compare retrieved accessions with those in the existing database to identify novel putative plasmid accessions that you may wish to include. Then, you could run the next stage of CuratePlasmids by providing the set of novel putative plasmids to the `--accessions` flag to determine plasmid accessions to be included in the existing database.



# Output files

The below table shows the outputs from running the complete pipeline (including retrieval of accessions from NCBI) with default settings.

File/Directory               | Description
---------------------------- | -------------------------------------------------------------------------------------------------
accessions.tsv               | putative plasmid accessions, following initial retrieval from NCBI
incompleteaccessions.tsv     | as above but not annotated as complete; these accessions will not be included
accessions_filtered.tsv      | accessions remaining, following filtering based on accession title text
accessions_filtered.fa       | sequences of accessions_filtered.tsv
accessions_excluded.tsv	     | accessions excluded, following filtering based on accession title text
plasmidfinder                | outputs from BLASTing accessions_filtered.fa against the PlasmidFinder replicon database
rmlst                        | outputs from BLASTing accessions_filtered.fa against the rMLST database
plasmids.fa		     | plasmid sequences
plasmids.tsv		     | plasmid accessions
rmlstrepaccessions.tsv	     | filtered accessions that have one or more replicon and rmlst loci detected (and are therefore not included in plasmids.fa)
rmlstonlyaccessions.tsv	     | filtered accessions that have one or more rmlst loci detected (and therefore are not included in plasmids.fa)

The .tsv output files contain the following columns: accession, topology (circular/linear), length, title, completeness (all are complete by default), rMLST loci (if non-plasmid)<br>
Accessions in rmlstrepaccessions.tsv are likely to be chromids. Accessions in rmlstonlyaccessions.tsv may be chromid sequences that don't contain a known plasmid replicon locus; alternatively, they may be chromosomal sequences mis-annotated as plasmids.


# Background and methods

For background information on curating plasmids see recent papers: [Orlek _et al._ (2017)](https://www.ncbi.nlm.nih.gov/pubmed/28286183) and [Galata _et al._ (2018)](https://academic.oup.com/nar/article/47/D1/D195/5149885). A brief outline of the steps of the complete pipeline (including retrieval from NCBI) is given below:

1. Putative complete plasmid accessions, along with accompanying information such as accession title are downloaded from NCBI nucleotide. To be considered a putative plasmid, the accession must be annotated as "plasmid" and "complete". Only Refseq accessions are retrieved. By default, source organism can be bacteria of any taxon, and accession creation date limits are not set.
2. The accessions are filtered using a regular expression search of the accession title text. For example, titles including the words "gene", "transposon", or "synthetic vector" would be excluded.
3. The filtered accession sequences are downloaded as a FASTA file.
4. The sequences are BLASTed again the PlasmidFinder replicon database and the rMLST database.
5. If a sequence contains no rMLST loci then it is considered a plasmid and included in the plasmids.fa output file. Accessions with rMLST loci detected are recorded.

# FAQ

* **How does CuratePlasmids differ from previously published plasmid curation methods?**
I previously published a similar method for plasmid curation (Orlek _et al._ (2017)](https://www.ncbi.nlm.nih.gov/pubmed/28286183), but compared with CuratePlasmids, the methods in the paper differ in several key ways:
* I used MLST rather than rMLST to filter chromosomal accessions. MLST loci are more limited as a chromosomal marker since MLST schemes cover fewer taxa.
* I included both Refseq and Genbank accessions, deduplicating identical sequences. In contrast, CuratePlasmids retrieves only Refseq accessions (and does not deduplicate based on sequence identity). Refseq represents a higher quality dataset (e.g. all Refseq accessions undergo the same annotation pipeline). Refseq is already a non-redundant database.
* Accessions were not required to be annotated as complete in the paper, but for CuratePlasmids, only complete plasmids are included.<br>
PLSDB follows similar methods to Orlek _et al_ 2017 but uses rMLST. However, unlike CuratePlasmids, accessions with up to 5 rMLST loci are included in PLSDB, whereas an accessions with one or more rMLST loci is not considered a plasmid when using CuratePlasmids - however, such accessions are recorded; they may be misannotated chromosomal sequence or chromids.
* **Does CuratePlasmids include linear plasmids?**
Yes, like the other plasmid databases described above, CuratePlasmids includes linear plasmids. However, the linear annotation should be treated cautiously since a 'linear' plasmid could actually be a circular plasmid that just couldn't be circularised after assembly. 


# Acknowledgements

I am grateful to Dr Keith Jolley informing me about programmatic access to the rMLST database and for pointing me towards [ConFindr](https://olc-bioinformatics.github.io/ConFindr/) software, which implements programmatic access. The `database_setup.py` executable used in CuratePlasmids is based on a script from ConFindr.


# License

[MIT License](https://en.wikipedia.org/wiki/MIT_License)
