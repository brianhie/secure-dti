import gzip
import numpy as np
import random
import sys

np.random.seed(10)
random.seed(10)

BINARY_SCORE=True

def load_chems(fname):
    chems = {}
    with open(fname) as f:
        for line in f:
            fields = line.rstrip().split()
            chem_id = fields[0]
            ecfp = fields[2]# + fields[3]
            chems[chem_id] = ecfp
    return chems

def load_prots(fname):
    prots = {}
    with open(fname) as f:
        for line in f:
            fields = line.rstrip().split()
            prot_id = fields[0]
            bitvec = fields[1]
            prots[prot_id] = bitvec
    return prots

def load_interactions(fname, chems, prots):
    if fname.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    interactions = []
    with opener(fname) as f:
        for x, line in enumerate(f):
            if interact_fname.endswith('.gz'):
                line = line.decode('utf-8')
            fields = line.rstrip().split()
            chem_id = fields[0]
            if not chem_id in chems:
                continue
            prot_id = fields[1]
            if prot_id.startswith('9606.'):
                prot_id = prot_id[len('9606.'):]
            if not prot_id in prots:
                continue
            if BINARY_SCORE:
                score = 1
            else:
                score = float(fields[2]) / 1000.
            interactions.append((chem_id, prot_id, score))
        
    return interactions

if __name__ == '__main__':
    interact_fname = sys.argv[1]
    chem_fname = sys.argv[2]
    prot_fname = sys.argv[3]
    dir_name = sys.argv[4].rstrip('/')
    
    divchem = False
    if len(sys.argv) >= 6:
        if sys.argv[5].lower() == 'divchem':
            divchem = True
            
    chems = load_chems(chem_fname)
    prots = load_prots(prot_fname)

    interactions = load_interactions(interact_fname, chems, prots)
    shuffled = interactions[:]
    int_set = set([ (c, p) for c, p, _ in interactions ])

    X_file = open(dir_name + '/X.txt', 'w')
    y_file = open(dir_name + '/y.txt', 'w')

    n = len(interactions)

    for idx, interaction in enumerate(interactions):
        # Write positive example.
        chem_id, prot_id, score = interaction
        feature = chems[chem_id] + prots[prot_id]
        X_file.write(feature + '\n')
        y_file.write(str(score) + '\n')

        # Write negative example.
        n_collisions = 0
        while True:
            # Fisher-Yates shuffle
            shuf_idx = np.random.randint(idx, n)
            (c1, p1, _) = shuffled[idx]
            (c2, p2, _) = shuffled[shuf_idx]
            # Do an edge swap.
            if divchem:
                p1, p2 = p2, p1
            else:
                c1, c2 = c2, c1
            if idx == n-1:
                break
            # Ensure that negative examples really are negative.
            if not (c2, p2) in int_set:
                break
            if n_collisions > 20:
                break
            n_collisions += 1
            
        if n_collisions > 20:
            sys.stderr.write('WARNING: Exceeded collision threshold.\n')
            continue

        shuffled[shuf_idx] = (c1, p1, 1)
        feature = chems[c2] + prots[p2]
        X_file.write(feature + '\n')
        y_file.write('0\n')

    X_file.close()
    y_file.close()
