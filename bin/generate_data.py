import gzip
import numpy as np
import random
import sys

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

    if interact_fname.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    X_file = open(dir_name + '/X.txt', 'w')
    y_file = open(dir_name + '/y.txt', 'w')

    with opener(interact_fname) as f:
        chem_idx = np.random.choice(range(len(chems.keys())), size=1600000)
        prot_idx = np.random.choice(range(len(prots.keys())), size=1600000)
        chem_l = chems.keys()
        prot_l = prots.keys()
        
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
            score = 1#float(fields[2]) / 1000.

            feature = chems[chem_id] + prots[prot_id]
            
            X_file.write(feature + '\n')
            y_file.write(str(score) + '\n')

            if divchem:
                feature = chems[chem_id] + prots[prot_l[prot_idx[x]]]
            else:
                feature = chems[chem_l[chem_idx[x]]] + prots[prot_l[prot_idx[x]]]
            X_file.write(feature + '\n')
            y_file.write('0\n')
            
    X_file.close()
    y_file.close()
