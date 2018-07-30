import datetime
import matplotlib.pyplot as plt
import numpy as np
import random
from sklearn import metrics
from string import ascii_lowercase
import sys

N_HIDDEN = 1
LOSS = 'hinge'

def report_scores(X, y, W, b, act):
    y_true = []
    y_pred = []
    y_score = []
    
    for l in range(N_HIDDEN):
        if l == 0:
            act[l] = np.maximum(0, np.dot(X, W[l]) + b[l])
        else:
            act[l] = np.maximum(0, np.dot(act[l-1], W[l]) + b[l])

    if N_HIDDEN == 0:
        scores = np.dot(X, W[-1]) + b[-1]
    else:
        scores = np.dot(act[-1], W[-1]) + b[-1]

    predicted_class = np.zeros(scores.shape)
    if LOSS == 'hinge':
        predicted_class[scores > 0] = 1
        predicted_class[scores <= 0] = -1
        y = 2 * y - 1
    else:
        predicted_class[scores >= 0.5] = 1

    sys.stdout.write(str(datetime.datetime.now()) + ' | ')
    print('Batch accuracy: {}'
          .format(metrics.accuracy_score(
              y, predicted_class
          ))
    )
    
    y_true.extend(list(y))
    y_pred.extend(list(predicted_class))
    y_score.extend(list(scores))
    
    # Output aggregated scores.
    try:
        sys.stdout.write(str(datetime.datetime.now()) + ' | ')
        print('Accuracy: {0:.2f}'.format(
            metrics.accuracy_score(y_true, y_pred))
        )
        sys.stdout.write(str(datetime.datetime.now()) + ' | ')
        print('F1: {0:.2f}'.format(
            metrics.f1_score(y_true, y_pred))
        )
        sys.stdout.write(str(datetime.datetime.now()) + ' | ')
        print('Precision: {0:.2f}'.format(
            metrics.precision_score(y_true, y_pred))
        )
        sys.stdout.write(str(datetime.datetime.now()) + ' | ')
        print('Recall: {0:.2f}'.format(
            metrics.recall_score(y_true, y_pred))
        )
        sys.stdout.write(str(datetime.datetime.now()) + ' | ')
        print('ROC AUC: {0:.2f}'.format(
            metrics.roc_auc_score(y_true, y_score))
        )
        sys.stdout.write(str(datetime.datetime.now()) + ' | ')
        print('Avg. precision: {0:.2f}'.format(
            metrics.average_precision_score(y_true, y_score))
        )
    except Exception as e:
        sys.stderr.write(str(e))
        sys.stderr.write('\n')
        
    return y_true, y_pred, y_score


def load_model():
    W = [ [] for _ in range(N_HIDDEN + 1) ]
    for l in range(N_HIDDEN+1):
        W[l] = np.loadtxt('mpc/cache/test_P1_W{}_final.bin'.format(l))

     # Initialize bias vector with zeros.
    b = [ []  for _ in range(N_HIDDEN + 1) ]
    for l in range(N_HIDDEN+1):
        b[l] = np.loadtxt('mpc/cache/test_P1_b{}_final.bin'.format(l))

    # Initialize activations.
    act = [ [] for _ in range(N_HIDDEN) ]

    return W, b, act
    
if __name__ == '__main__':
    W, b, act = load_model()

    X_train = np.genfromtxt('data/batch_pw/Xtrain',
                            delimiter=1, dtype='float')
    y_train = np.genfromtxt('data/batch_pw/ytrain',
                            delimiter=1, dtype='float')
    
    print('Training accuracy:')
    report_scores(X_train, y_train, W, b, act)

    X_test = np.genfromtxt('data/batch_pw/Xtest',
                            delimiter=1, dtype='float')
    y_test = np.genfromtxt('data/batch_pw/ytest',
                            delimiter=1, dtype='float')
    
    print('Testing accuracy:')
    report_scores(X_test, y_test, W, b, act)
