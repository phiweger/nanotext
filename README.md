## nanotext

Name inspired by [`fastText`](https://fasttext.cc/).

Run tests: 


```bash
cd .../nanotext/
pytest  # or python setup.py test
```


```
TODO: turn embedding in gensim fmt into something similar to glove (bin?)
would then be easier to train/ unfreeze and then save again, to be loaded w/ some model eg for PUL prediction

https://github.com/plasticityai/magnitude

nanotext compute

takes annotation (load domains) and our model and computes vector

nanotext train corpus model

nanotext compare ...

nanotext search ...

nanotext taxonomy (calculate a gtdb based taxonomy and get closest functional genome and use that or distance)

like sourmash really



check sourmash publication

nanotext predict model=medium
nanotext predict model=pul
```