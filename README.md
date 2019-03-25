# nanotext

This library enables the use of embedding vectors generated from a large corpus of protein domains to search for _similar_ genomes, where similar is the cosine similarity between one genome's vector and another's. Think about protein domains as words, genomes as documents, and search as a form of document retrieval based on the notion of _topic_. You can read more about this work in out [bioRxiv preprint](https://www.biorxiv.org/content/early/2019/01/18/524280). Name inspired by [`fastText`](https://fasttext.cc/).

All releases of `nanotext` starting from r89 are pegged against the corresponding release of the [Genome Taxonomy Database (GTDB)](http://gtdb.ecogenomic.org/).


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


## Data

We provide all data through the [Open Science Framework (OSF)](https://osf.io) because [reasons](http://ivory.idyll.org/blog/2017-osf-for-files.html). The project ID is [`pjf7m`](https://osf.io/pjf7m/) and you can either download the files manually or use the OSF client:


```bash
osf -p pjf7m fetch <filename>
```


## Quick start

Head over to the [tutorial](https://github.com/phiweger/nanotext/blob/master/tutorial/tara.ipynb) for a quick walkthrough.


## Training

If you're keen, you can [download the corpus](https://osf.io/pjf7m/) (6 GB unpacked) and train the embedding yourself (couple of hours, about 10 GB of RAM). Note that training is stochastic, and that the exact similarity values will differ slightly from the numbers below. The current release of `nanotext` (r89), incorporates about 145 thousand genomes with one billion domains from the [Genome Taxonomy Database (GTDB)](http://gtdb.ecogenomic.org/). The associated corpus and models are available from OSF. We provide a [training and evaluation workflow](https://github.com/phiweger/nanotext/tree/master/nanotext/workflows/model_training).


## Search similar genomes and infer taxonomy

As an example, let's pass in a genome annotation of a _Prochlorococcus_ genome assembly based on data from the _Tara Ocean Expedition_ by [Delmont, T. O. et al. Nat Microbiol 3, 804–813 (2018)](https://www.nature.com/articles/s41564-018-0176-9).

If you don't have annotations for your favourite microbe already, we provided a small `snakemake` workflow [here](https://github.com/phiweger/nanotext/tree/master/nanotext/workflows/annotation_pfamscan).

We'll use the pretrained embedding to search for _functionally_ similar genomes, where function is analogous to the topic of a document.


```bash
osf -p pjf7m fetch models.zip && unzip models.zip
osf -p pjf7m fetch tara.zip && unzip tara.zip
osf -p pjf7m fetch metadata_GTDB_r89.db

nanotext search --models models --topn 30 --mode core \
  --annotation tara/TARA_ION_MAG_00012_pfam.tsv \
  --out TARA_ION_MAG_00012.most_similar.tsv

head -n3 TARA_ION_MAG_00012.most_similar.tsv
# GCF_000158595.1 0.9534
# GCF_000760055.1 0.9481
# GCA_003281365.1 0.9462
```


Sometimes it's interesting to know which taxa these similar genomes are from, e.g. when trying to identify MAGs.


```bash
nanotext search --models models 
  --annotation tara/TARA_ION_MAG_00012_pfam.tsv \
  --taxonomy metadata_GTDB_r89.db
# will produce taxonomy.tsv

head -n2 taxonomy.tsv
# name  domain  phylum  class order family  genus species
# GCF_000158595.1 Bacteria  Cyanobacteriota Cyanobacteriia  Synechococcales_A CyanobiaceaeProchlorococcus_A Prochlorococcus_A sp5
```


## Prediction

You can use genome vectors as direct input to machine learning algorithms. We provide a prove of principle, predicting culture media from a genome's protein domain annotation alone. This model was trained on the GTDB r83 and will NOT be ported to future releases. However, you can query your own MAGs nonetheless. Note also that for historical reasons, the annotation will have to be using `hmmer`, for which we provided a `snakemake` workflow [here](https://github.com/phiweger/nanotext/tree/master/nanotext/workflows/annotation_hmmer).


```bash
osf -p pjf7m fetch culture.zip && unzip culture.zip
nanotext predict \
    --embedding culture/embedding.genomes.model \
    --genome culture/TARA_ION_MAG_00012.domtbl.tsv \
    --model culture/media_prediction.h5 \
    --db culture/embedding.media.json \
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

