## Model training and evaluation

We assume that for each set of hyperparameters you want to train, a corresponding config file is present in the workflow's `config/nanotext/` folder.

Two evaluation tasks are currently implemented:

- Semantiv odd-one out (SOMO) task
- Ecotype task

Please refer to the preprint for details.


```bash
# First, download and unzip corpus from OSF, e.g. using osfclient:
mkdir results && cd results
osf -p pjf7m fetch corpus_r89.txt.zip && unzip corpus_89.txt.zip

# Train the embedding
# (1) cd into the workflow directory
# (2) on a Mac, to avoid AppNap (training stopping halfway through) run
defaults write org.python.python NSAppSleepDisabled -bool YES
# (3) set the paths in config/snakemake.json
# (4) train already
snakemake -p --configfile config/snakemake.json --cores 8
```