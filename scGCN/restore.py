from __future__ import division
from __future__ import print_function
import os
import sys
import time
import numpy as np
import pickle as pkl
import tensorflow as tf
from utils import *
from tensorflow.python.saved_model import tag_constants
from models import *
#' del_all_flags(FLAGS)

# Set random seed
seed = 123
np.random.seed(seed)
tf.compat.v1.set_random_seed(seed)
tf.set_random_seed(seed)

# Settings
flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_string('dataset', 'Data_gathered', 'data dir')
flags.DEFINE_string('model', 'scGCN','Model string.') 
                     
# Load data
adj, features, labels_binary_train, labels_binary_val, labels_binary_test, train_mask, pred_mask, val_mask, test_mask, true_label = load_data(
    FLAGS.dataset)

#' load results 
import pickle as pkl
PIK = './results/results.dat'
with open(PIK, 'rb') as f:
    objects = pkl.load(f)

train_accuracy, test_accuracy, val_accuracy, train_loss, test_loss, val_loss, activation_output, predict_output, acc_pred = tuple(
    objects)

#' accuracy on all masks 
labels_binary_all = new_label

sess = tf.Session()
ab = sess.run(tf.nn.softmax(predict_output))
all_prediction = sess.run(
    tf.equal(sess.run(tf.argmax(ab, 1)),
             sess.run(tf.argmax(labels_binary_all, 1))))

all_mask = np.array([True] * len(train_mask))
labels_binary_all = np.concatenate(
    (labels_binary_train[train_mask], labels_binary_val[val_mask],
     labels_binary_test[test_mask]))

y_pred = sess.run(tf.argmax(ab, 1))[pred_mask]
y_real = sess.run(tf.argmax(labels_binary_all, 1))[pred_mask]

from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score
print("Prediction accuracy :{}".format(acc_pred))
print("F1 score :{}".format(f1_score(y_real, y_pred, average='weighted')))
print("Precision score :{}".format(
    precision_score(y_real, y_pred, average='weighted')))
print("Recall score :{}".format(
    recall_score(y_real, y_pred, average='weighted')))
