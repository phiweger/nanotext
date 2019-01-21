# Workflows

To test the annotation workflow, we provided a test MAG in `genomes/`. To run, first modify the `/path/to/pfam_v32/Pfam-A.hmm` in `config.json`. Then just execute the following with this folder as the working directory:


```bash
# download and set up HMM database
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm

# set config path to wherever the HMM database is, then:
snakemake --configfile config.json --cores=8  # make this suit your machine
```


You should see a folder `results/` being created with the annotation in `data/` and a `log/` if anything goes wrong. Note that when you run this workflow on your own genomes and their file extension is not `.fa`, you need to modify the Snakefile accordingly.