# nanotext

This library enables the use of embedding vectors generated from a large corpus of protein domains to search for _similar_ genomes, where similar is the cosine similarity between one genome's vector and another's. Think about protein domains as words, genomes as documents, and search as a form of document retrieval based on the notion of _topic_. 

You can read more about this work in out [bioRxiv preprint](https://www.biorxiv.org/content/early/2019/01/18/524280). Name inspired by [`fastText`](https://fasttext.cc/).

Install and run tests: 


```bash
git clone https://github.com/phiweger/nanotext
cd nanotext
pip install -e .
pytest  # or python setup.py test
```


As an example, pass in a genome annotation of a _Prochlorococcus_ genome assembly based on data from the _Tara Ocean Expedition_ by [Delmont, T. O. et al. Nat Microbiol 3, 804–813 (2018)](https://www.nature.com/articles/s41564-018-0176-9).

If you don't have annotations already, we provided a small `snakemake` workflow [here](https://github.com/phiweger/nanotext/tree/master/nanotext/workflows/annotation_tara).

We'll use the pretrained embedding to search for _functionally_ similar genomes, where function is analogous to the topic of a document. If you're keen, you can [download the corpus](https://osf.io/pjf7m/) (280 MB) and train the embedding yourself (about 3 hours and 10 GB RAM on a recent MacBook Pro). Note that training is stochastic, and that the exact similarity values will differ slightly from the numbers below.


```bash
# first download and unzip corpus; then to train the embedding:
nanotext train --corpus data/corpus.txt --out data/embedding.genomes.model

# or use pretrained model
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

nanotext predict \\
    --embedding data/embedding.genomes.model \\
    --genome data/TARA_ION_MAG_00012.domtbl.tsv \\
    --model data/media_prediction.h5 \\
    --db data/embedding.media.json \\
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
