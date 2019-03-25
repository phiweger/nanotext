# nanotext

This library enables the use of embedding vectors generated from a large corpus of protein domains to search for _similar_ genomes, where similar is the cosine similarity between one genome's vector and another's. Think about protein domains as words, genomes as documents, and search as a form of document retrieval based on the notion of _topic_. 

You can read more about this work in out [bioRxiv preprint](https://www.biorxiv.org/content/early/2019/01/18/524280). Name inspired by [`fastText`](https://fasttext.cc/).

Starting from "r89", all future releases of `nanotext` are pegged against the corresponding release of the [Genome Taxonomy Database (GTDB)](http://gtdb.ecogenomic.org/).


## Installation and tests


```bash
# Create a new conda environment to experiment in
conda create -y -n myenv && conda activate myenv
conda install faiss-cpu -c pytorch  # for fast vector search

# Install the nanotext library
git clone https://github.com/phiweger/nanotext
cd nanotext
pip install -e .
pytest  # or python setup.py test

# For the tutorial, we need some notebook and visualisation libraries
conda install jupyter
conda install -c conda-forge altair vega_datasets notebook vega

# To download from the project data on OSF using the command line
conda install -c conda-forge osfclient
```


## Quick start

Head over to the [tutorial](https://github.com/phiweger/nanotext/blob/master/tutorial/tara.ipynb) for a quick walkthrough.


## Training

Update: We are working on a new release of `nanotext` (r89), incorporating about 145k genomes from the [Genome Taxonomy Database (GTDB)](http://gtdb.ecogenomic.org/). The associated corpus and models are available from OSF, project ID `pjf7m` (see the [tutorial](https://github.com/phiweger/nanotext/blob/master/tutorial/tara.ipynb)).

We provide a training and evaluation workflow. It assumes that for each set of hyperparameters you want to train, a corresponding config file is present is the workflow's `config/nanotext/` folder.


```bash
# First, download and unzip corpus from OSF, e.g. using osfclient:
mkdir results && cd results
osf -p pjf7m fetch corpus_r89.txt.zip corpus_r89.txt.zip
unzip corpus_89.txt.zip

# Train the embedding
# (1) cd into the workflow directory
# (2) on a Mac, to avoid AppNap (training stopping halfway through) run
defaults write org.python.python NSAppSleepDisabled -bool YES
# (3) train already
snakemake -p --configfile config/snakemake.json --cores 8
```


## Search similar genomes and infer taxonomy

As an example, pass in a genome annotation of a _Prochlorococcus_ genome assembly based on data from the _Tara Ocean Expedition_ by [Delmont, T. O. et al. Nat Microbiol 3, 804–813 (2018)](https://www.nature.com/articles/s41564-018-0176-9).

If you don't have annotations already, we provided a small `snakemake` workflow [here](https://github.com/phiweger/nanotext/tree/master/nanotext/workflows/annotation_hmmer).

We'll use the pretrained embedding to search for _functionally_ similar genomes, where function is analogous to the topic of a document. If you're keen, you can [download the corpus](https://osf.io/pjf7m/) (280 MB) and train the embedding yourself (about 3 hours and 10 GB RAM on a recent MacBook Pro). Note that training is stochastic, and that the exact similarity values will differ slightly from the numbers below.


## Prediction

...


---

DEPRECATED, will be replaced soon ...


```bash
# Alternatively, use the pretrained model provided in this repo:
nanotext search \
    --embedding data/embedding.genomes.model \
    --genome data/TARA_ION_MAG_00012.domtbl.tsv \
    --topn 3 --out -
# Loading model ...
# Inferring genome vector ...
# GCA_000634215.1 0.9344
# GCF_000759935.1 0.9282
# GCF_000759855.1 0.9276
# Done.

nanotext predict \
    --embedding data/embedding.genomes.model \
    --genome data/TARA_ION_MAG_00012.domtbl.tsv \
    --model data/media_prediction.h5 \
    --db data/embedding.media.json \
    --topn 3 --out -
# Using TensorFlow backend.
# Loading embedding model for genomes ...
# Inferring genome vector ...
# Loading embedding model for culture media ...
# Loading predictive model ...
# 137 0.9915
# 918 0.987
# 69  0.986
# Done.
```


The result is a list of media IDs and their cosine similarity to the prediction. Now you can check out the associated media ingredients from the [DSMZ list of recommended media for microorganisms](https://www.dsmz.de/catalogues/catalogue-microorganisms/culture-technology/list-of-media-for-microorganisms.html).


## Taxonomy

Given a contig from a metagenome assembly, you can query its likely taxonomy like so:


```bash
nanotext taxonomy \
    --embedding nanotext_r89.model --taxonomy bac_taxonomy_r86.tsv \
    --query JFOD01_pfam.tsv --fmt pfamscan --topn 10 -o results.json

cat results.json | jq ".majority,.ratio"
```


This is what the result looks like (the file contains raw taxonomy and cosine similarity values for the top hits, too, in case this is of value):



```json
{
  "domain": "Bacteria",
  "phylum": "Cyanobacteriota",
  "class": "Cyanobacteriia",
  "order": "Synechococcales_A",
  "family": "Cyanobiaceae",
  "genus": "Prochlorococcus_A",
  "species": "Prochlorococcus_A sp1"
}
{
  "domain": 1,
  "phylum": 1,
  "class": 1,
  "order": 1,
  "family": 1,
  "genus": 0.83,
  "species": 0.5
}
```
