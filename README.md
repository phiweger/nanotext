## nanotext

Library that uses embedding vectors generated from a large corpus of protein domains to search for _similar_ genomes, where similar is the cosine similarity between one genome's vector and another's. Think about protein domains as words, genomes as documents, and search as a form of document retrieval based on the notion of _topic_.

Name inspired by [`fastText`](https://fasttext.cc/).

Install and run tests: 


```bash
git clone https://github.com/phiweger/nanotext
cd nanotext
pytest  # or python setup.py test
```


As an example, pass in a genome annotation of a _Prochlorococcus_ genome assembly based on data from the _Tara Ocean Expedition_ by 

Delmont, T. O. et al. Nat Microbiol 3, 804–813 (2018)

We'll use the pretrained embedding to search for _functionally_ similar genomes, where function is analogous to the topic of a document.


```bash
# download corpus from OSF.io and train embedding
nanotext train --corpus data/corpus.txt --out data/embedding.genomes.model

# or use pretrained model
nanotext search \
    --embedding data/embedding.genomes.model \
    --genome data/TARA_ION_MAG_00012.domtbl.tsv \
    --topn 3 --out -
# Loading model ...
# Inferring genome vector ...
# GCA_000634215.1 0.93
# GCF_000759935.1 0.93
# GCF_000759855.1 0.93
# Done.
```


More to follow, watch this space.