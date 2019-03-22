from pathlib import Path

import faiss
import numpy as np
from sklearn.preprocessing import normalize

from nanotext.utils import subtract_mean, index_model, infer_genome_vector
from nanotext.io import eprint, load_embedding


class GenomeModel():
    '''
    No norm means distance is assumed to be Euclidean.

    Usage:

    # load
    from nanotext.classes import GenomeModel
    fp_models = 'path/to/gensim/models'
    ensemble = GenomeModel(fp_models, mode='core', norm='l2')

    # example query
    fp_query = 'path/to/query.pfam.tsv'
    v = ensemble.infer(fp_query, fmt='pfamscan')
    ensemble.search(v, topn=5, min_dist=0.5)

    # lookup
    from nanotext.utils import cosine
    cosine(ensemble['GCA_003529605.1'], ensemble['GCA_002433265.1'])
    '''
    def __init__(self, fp, mode='ensemble', norm=None, names=None):
        
        self.mode = mode
        if mode == 'ensemble':
            nn = ['22', '45', '93']
        elif mode == 'core':
            nn = ['93']
        elif mode == 'accessory':
            nn = ['22']
        else:
            raise ValueError(
                'More not implemented (try "ensemble", "core" or "accessory")')
        
        self.models = []
        for n in nn:
            p = Path(fp) / f'{n}/nanotext_r89.model'
            model = load_embedding(str(p))
            self.models.append(model)
    
        self.norm = norm
        eprint('Subtracting mean from model(s) ...')
        self.nomean = self._demean(self.models)

        if not names:
            self.names = self.models[0].docvecs.index2entity

        self.dim = len(self.models[0].docvecs[0])
        
        eprint('Indexing model(s) ...')
        if norm:
            eprint(f'{self.norm} norm will be applied to vectors')
        found, m, self.index = index_model(
            self.names, [i for i in self.nomean], self.norm)
        self.embedding = dict(zip(found, m))
    
        self.means = []
        for model in self.models:
            mu = np.mean(
                [model.docvecs[i] for i in range(len(model.docvecs))], axis=0)
            self.means.append(mu)

        self.warn_on_ensemble_inference = False


    def _demean(self, models):
        for i in models:
            yield subtract_mean(i)


    def infer(self, fp, steps=1000, fmt='pfamscan', truncate_by=0):
        if (self.mode == 'ensemble') and (not self.warn_on_ensemble_inference):
            eprint('''Warning: Inference w/ a model ensemble will work well if you don't combine the resulting vectors w/ the indexed ones. This is because small variations in the inference will magnify in model ensembles to offset the inferred and indexed vectors by more than they actually differ.
                ''')
            self.warn_on_ensemble_inference = True  # print only once
        bag = []

        for model, mu in zip(self.models, self.means):
            v = infer_genome_vector(
                fp, model, steps=steps, fmt=fmt, truncate_by=truncate_by)
            v_ = v-mu

            bag.append(v_)
        
        ve = np.mean(bag, axis=0)  # ensemble vector
        ve = np.array([ve], dtype='float32')  # cast for norm and index search
        
        if self.norm == 'l2':
            # eprint('L2 normalization ...')
            ve = normalize(ve, norm=self.norm, axis=1)#.reshape(self.dim)
            # w/o reshape, dim is (1, dim), not (dim,) like the model's vecs;
            # this causes problems when we want to combine trained and inferred
            # vecs
            # on the diff btw/ (100,) and (100, 1) see
            # stackoverflow.com/questions/22053050

        return ve
        # so we can directly use it w/ search 


    def search(self, query, topn=3, min_dist=None):
        D, I = self.index.search(query, topn)
        if not min_dist:
            return [(self.names[i], j) for i, j in zip(I[0], D[0])]
        else:
            return [(self.names[i], j) for i, j in zip(I[0], D[0]) \
                    if j>min_dist]


    def subset(self, names):
        '''
        Given a list of names, return a dict of name: vector
        '''
        d = self.embedding.copy()
        embedding = {}
        notfound = 0
        for i in names:
            try:
                v = d[i]
                embedding[i] = v
            except KeyError:
                notfound += 1
        eprint(f'{notfound} records not found')
        return embedding

    
    def __getitem__(self, key):
        return self.embedding[key].squeeze()
