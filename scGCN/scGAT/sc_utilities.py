import pickle as pkl
import scipy.sparse
import numpy as np
import pandas as pd
from scipy import sparse as sp
import networkx as nx
from gcn.utils import *
from sc_data import *
from collections import defaultdict
from scipy.stats import uniform
#' -------- convert graph to specific format -----------


def get_value(diction, specific):
    for key, val in diction.items():
        if val == specific:
            return (key)


def graph(matrix):
    adj = defaultdict(list)  # default value of int is 0
    for i, row in enumerate(matrix):
        for j, adjacent in enumerate(row):
            if adjacent:
                adj[i].append(j)
        if adj[i].__len__ == 0:
            adj[i] = []
    return adj


def sample_mask(idx, l):
    """Create mask."""
    mask = np.zeros(l)
    mask[idx] = 1
    return np.array(mask, dtype=np.bool)


# convert nested lists to a flat list
output = []


def removNestings(l):
    for i in l:
        if type(i) == list:
            removNestings(i)
        else:
            output.append(i)
    return (output)


def load_data(datadir):
    input_data(datadir)
    PIK = "{}/datasets.dat".format(datadir)
    with open(PIK, "rb") as f:
        objects = pkl.load(f)

    data_train1, data_test1, data_val1, label_train1, label_test1, label_val1, lab_data2, lab_label2, types = tuple(
        objects)

    train2 = pd.concat([data_train1, lab_data2])
    lab_train2 = pd.concat([label_train1, lab_label2])

    datas_train = np.array(train2)
    datas_test = np.array(data_test1)
    datas_val = np.array(data_val1)

    labels_train = np.array(lab_train2).flatten()
    labels_test = np.array(label_test1).flatten()
    labels_val = np.array(label_val1).flatten()

    #' convert pandas data frame to csr_matrix format
    datas_tr = scipy.sparse.csr_matrix(datas_train.astype('Float64'))
    datas_va = scipy.sparse.csr_matrix(datas_val.astype('Float64'))
    datas_te = scipy.sparse.csr_matrix(datas_test.astype('Float64'))

    #' 3) set the unlabeled data in training set

    #' @param N; the number of labeled samples in training set
    M = len(data_train1)

    #' 4) get the feature object by combining training, test, valiation sets

    features = sp.vstack((sp.vstack((datas_tr, datas_va)), datas_te)).tolil()
    #'features = preprocess_features(features)

    #' 5) Given cell type, generate three sets of labels with the same dimension
    labels_tr = labels_train.flatten()
    labels_va = labels_val.flatten()
    labels_te = labels_test.flatten()

    labels = np.concatenate(
        [np.concatenate([labels_tr, labels_va]), labels_te])
    Labels = pd.DataFrame(labels)

    true_label = Labels
    #' convert list to binary matrix
    uniq = np.unique(Labels.values)

    rename = {}

    for line in range(0, len(types)):
        key = types[line]
        rename[key] = int(line)

    Label1 = Labels.replace(rename)
    indices = np.array(Label1.values, dtype='int').tolist()

    indice = [item for sublist in indices for item in sublist]

    #' convert list to binary matrix
    indptr = range(len(indice) + 1)
    dat = np.ones(len(indice))
    binary_label = scipy.sparse.csr_matrix((dat, indice, indptr))

    #' new label with binary values
    new_label = np.array(binary_label.todense())
    idx_train = range(M)
    idx_pred = range(M, len(labels_tr))
    idx_val = range(len(labels_tr), len(labels_tr) + len(labels_va))
    idx_test = range(
        len(labels_tr) + len(labels_va),
        len(labels_tr) + len(labels_va) + len(labels_te))

    train_mask = sample_mask(idx_train, new_label.shape[0])
    pred_mask = sample_mask(idx_pred, new_label.shape[0])
    val_mask = sample_mask(idx_val, new_label.shape[0])
    test_mask = sample_mask(idx_test, new_label.shape[0])

    labels_binary_train = np.zeros(new_label.shape)
    labels_binary_val = np.zeros(new_label.shape)
    labels_binary_test = np.zeros(new_label.shape)
    labels_binary_pred = np.zeros(new_label.shape)
    labels_binary_train[train_mask, :] = new_label[train_mask, :]
    labels_binary_val[val_mask, :] = new_label[val_mask, :]
    labels_binary_test[test_mask, :] = new_label[test_mask, :]
    labels_binary_pred[pred_mask, :] = new_label[pred_mask, :]

    #' ----- use seurat output to construct matrix ---------

    id_graph1 = pd.read_csv('{}/inter_graph.csv'.format(datadir),
                            index_col=0,
                            sep=',')
    id_graph2 = pd.read_csv('{}/intra_graph.csv'.format(datadir),
                            sep=',',
                            index_col=0)

    #' --- map index ----
    fake1 = np.array([-1] * len(lab_data2.index))
    index1 = np.concatenate((data_train1.index, fake1, data_val1.index,
                             data_test1.index)).flatten()
    #' (feature_data.index==index1).all()
    fake2 = np.array([-1] * len(data_train1))
    fake3 = np.array([-1] * (len(data_val1) + len(data_test1)))
    find1 = np.concatenate((fake2, np.array(lab_data2.index), fake3)).flatten()

    #' ---------------------------------------------
    #'  graph 2
    #' ---------------------------------------------
    id_grp1 = np.array([
        np.concatenate((np.where(find1 == id_graph2.iloc[i, 1])[0],
                        np.where(find1 == id_graph2.iloc[i, 0])[0]))
        for i in range(len(id_graph2))
    ])

    id_grp2 = np.array([
        np.concatenate((np.where(find1 == id_graph2.iloc[i, 0])[0],
                        np.where(find1 == id_graph2.iloc[i, 1])[0]))
        for i in range(len(id_graph2))
    ])

    #' ---------------------------------------------
    #'  inter-graph
    #' ---------------------------------------------
    id_gp1 = np.array([
        np.concatenate((np.where(find1 == id_graph1.iloc[i, 1])[0],
                        np.where(index1 == id_graph1.iloc[i, 0])[0]))
        for i in range(len(id_graph1))
    ])

    id_gp2 = np.array([
        np.concatenate((np.where(index1 == id_graph1.iloc[i, 0])[0],
                        np.where(find1 == id_graph1.iloc[i, 1])[0]))
        for i in range(len(id_graph1))
    ])

    matrix = np.identity(len(labels))
    matrix[tuple(id_grp1.T)] = 1
    matrix[tuple(id_grp2.T)] = 1
    matrix[tuple(id_gp1.T)] = 1
    matrix[tuple(id_gp2.T)] = 1

    adj = graph(matrix)
    adj = nx.adjacency_matrix(nx.from_dict_of_lists(adj))

    print("assign input coordinatly....")
    return adj, features, labels_binary_train, labels_binary_val, labels_binary_test, labels_binary_pred, train_mask, pred_mask, val_mask, test_mask, new_label, true_label


#' in case FLAGS are defined twice
def del_all_flags(FLAGS):
    flags_dict = FLAGS._flags()
    keys_list = [keys for keys in flags_dict]
    for keys in keys_list:
        FLAGS.__delattr__(keys)
