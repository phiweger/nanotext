def chunks(l, n):
    '''
    Yield successive n-sized chunks from l (stackoverflow, 312443).

    a = [1, 2, 3, 4]
    list(chunks(a, 2))
    # [[1, 2], [3, 4]]

    Returns empty list if list empty.

    For overlapping chunks, see windows()
    '''
    for i in range(0, len(l), n):
        yield l[i:i + n]


def odd_one(fp_test_corpus, fp_model, n_not_odd=5):
    '''
    SOMO task.
    '''
    import random
    from nanotext.io import load_embedding

    model = load_embedding(fp_model)
    vocab = list(model.wv.vocab.keys())
    pos, neg = 0, 0

    with open(fp_test_corpus, 'r') as file:
        for line in file:
            genome, contig, domains = line.strip().split('\t')
            domains = domains.split(',')
            for seq in chunks(domains, n_not_odd):
                if len(seq) > 1:  # otherwise its a 50:50 coin flip
                    odd = random.choice(vocab)
                    guess = model.wv.doesnt_match(seq+[odd])
                    if guess == odd:
                        pos += 1
                    else:
                        neg += 1

    return round(pos/(pos+neg), 4)










'''
# load model and vocabulary
model = Doc2Vec.load('embedding.doc2vec.model')
vocab = list(model.wv.vocab.keys())  # 10879 domains


# load list of filepaths from Ensembl
fp = '/Volumes/container/nanotext/data/ensembl_annotation/ftp.ensemblgenomes.org/pub/bacteria/release-35/json'
annos = glob(f'{fp}/**/*.json', recursive=True)


print('Creating validation set from annotation files ...')
validation, cnt = [], 0

for p in tqdm(annos):  # p .. path
    anno = is_test(p)
    if anno:
        cnt += 1
        print(cnt)
        if cnt == 100:
            break

        for g in anno['genes']:
            try:
                pfam = g['Pfam']  # get an ORF's domains
                if len(pfam) > 1:                         # multidomain protein
                    pfam.append(random.choice(vocab))     # add the "odd one"
                    validation.append(pfam)               # last one is odd one
                else:
                    pass                 
            except KeyError:
                pass


count, correct = 0, 0
# wrong = []
for i in validation:
    count += 1
    odd = i[-1]  # we appended the random one to the last position
    guess = model.wv.doesnt_match(i)
    # How this thing works: stackoverflow.com/questions/45948533
    if guess == odd:
        correct += 1
        # print(f'\t\t\t\tcorrect: {i[-1]}, called: {guess}')
    else:
        # print(f'correct: {i[-1]}, called: {guess}')
        # wrong.append(odd)
        pass

# print(sorted(wrong))
print(count, correct)
print(f'The embedding has an accuracy of {round(correct/count, 4)}')
# The embedding has an accuracy of 0.9927
# 0.9926664103853091
len(test_samples)  # 873
'''


def semantics():
    '''
    king queen
    '''
    pass


def cluster():
    pass

