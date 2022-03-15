import pickle as pkl
import scipy.sparse
import numpy as np
import pandas as pd
from scipy import sparse as sp
import networkx as nx
from data import *
from collections import defaultdict
from scipy.stats import uniform
import tensorflow as tf

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
    # 给定性状和类型的用0填充的数组
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

def load_data(datadir,rgraph=True):
    input_data(datadir,Rgraph=rgraph)
    PIK = "{}/datasets.dat".format(datadir)
    with open(PIK, "rb") as f:
        objects = pkl.load(f)
    data_train1, data_test1, data_val1, label_train1, label_test1, label_val1, lab_data2, lab_label2, types = tuple(
        objects)
    # data_train1 data_test1 data_val1只有ref,label相对应
    # lab_data2与lab_label2是qur
    #此处原来的types仅仅是ref中的types
    train2 = pd.concat([data_train1, lab_data2])
    lab_train2 = pd.concat([label_train1, lab_label2])
    #train2包含了两个数据集的数据，包含了从Ref中随机得到的train以及所有的qur
    datas_train = np.array(train2)
    datas_test = np.array(data_test1)
    datas_val = np.array(data_val1)
    ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ####!!!!!!!!!!!!!!!!!!!!这里就是给未来的结果分配index,依次是ref_train/qur/ref_val/ref_test
    index_guide = np.concatenate(
        (label_train1.index, lab_label2.index * (-1) - 1, label_val1.index,
         label_test1.index))
    #numpy提供了numpy.concatenate((a1,a2,...), axis=0)函数。能够一次完成多个数组的拼接
    labels_train = np.array(lab_train2).flatten()
    labels_test = np.array(label_test1).flatten()
    labels_val = np.array(label_val1).flatten()
    #同上，labels_train包含了两个数据集
    #' convert pandas data frame to csr_matrix format
    datas_tr = scipy.sparse.csr_matrix(datas_train.astype('Float64'))
    datas_va = scipy.sparse.csr_matrix(datas_val.astype('Float64'))
    datas_te = scipy.sparse.csr_matrix(datas_test.astype('Float64'))

    #' 3) set the unlabeled data in training set
    #' @param N; the number of labeled samples in training set
    M = len(data_train1)

    #这里是有改动？
    #' 4) get the feature object by combining training, test, valiation sets
    features = sp.vstack((sp.vstack((datas_tr, datas_va)), datas_te)).tolil()
    features = preprocess_features(features)

    #' 5) Given cell type, generate three sets of labels with the same dimension
    labels_tr = labels_train.flatten()
    labels_va = labels_val.flatten()
    labels_te = labels_test.flatten()

    labels = np.concatenate(
        [np.concatenate([labels_tr, labels_va]), labels_te])
    Labels = pd.DataFrame(labels)

    true_label = Labels
    #这个变量中存储了所有的原本真实标签label
    #' convert list to binary matrix
    uniq = np.unique(Labels.values)

    rename = {}
    #这里把类型转换成了数字，而此处的类型只有ref的;原始types变量里面只存储了唯一的Ref的label
    for line in range(0, len(types)):
        key = types[line]
        rename[key] = int(line)
    #！！！！这里把类型数字进行了替换，所以会导致data2内新类型仍然是文字而报错
    Label1 = Labels.replace(rename)
    #Labels是原始的所有的真实标签label
    #因此Lanel1是所有转换为数字？的真实标签
    indices = np.array(Label1.values, dtype='int').tolist()

    indice = [item for sublist in indices for item in sublist]

    #' convert list to binary matrix
    indptr = range(len(indice) + 1)
    dat = np.ones(len(indice))
    binary_label = scipy.sparse.csr_matrix((dat, indice, indptr))
    #binary_label是所有的二进制标签
    #' new label with binary values
    new_label = np.array(binary_label.todense())
    #new_label是所有的二进制的值
    idx_train = range(M)
    #M = len(data_train1) data_train1是从ref中提取出来的train
    idx_pred = range(M, len(labels_tr))
    #这个区间相当于是所有的qur
    idx_val = range(len(labels_tr), len(labels_tr) + len(labels_va))
    #这个区间仅仅是从ref中提取出来的val
    idx_test = range(
        len(labels_tr) + len(labels_va),
        len(labels_tr) + len(labels_va) + len(labels_te))
    #这个区间仅仅是从ref中提取出来的test
    train_mask = sample_mask(idx_train, new_label.shape[0])
    pred_mask = sample_mask(idx_pred, new_label.shape[0])
    val_mask = sample_mask(idx_val, new_label.shape[0])
    test_mask = sample_mask(idx_test, new_label.shape[0])
    #给定形状和类型的用0填充的数组，也就是所有的标签做01处理
    labels_binary_train = np.zeros(new_label.shape)
    labels_binary_val = np.zeros(new_label.shape)
    labels_binary_test = np.zeros(new_label.shape)
    labels_binary_train[train_mask, :] = new_label[train_mask, :]
    labels_binary_val[val_mask, :] = new_label[val_mask, :]
    labels_binary_test[test_mask, :] = new_label[test_mask, :]

    #' ----- construct adjacent matrix ---------

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
    #这里是什么意思啊？师兄？
    #' ---------------------------------------------
    #'  intra-graph
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
    return adj, features, labels_binary_train, labels_binary_val, labels_binary_test, train_mask, pred_mask, val_mask, test_mask, new_label, true_label, index_guide, rename 


def preprocess_features(features):
    """Row-normalize feature matrix and convert to tuple representation"""
    rowsum = np.array(features.sum(1))
    r_inv = np.power(rowsum, -1).flatten()
    r_inv[np.isinf(r_inv)] = 0.
    r_mat_inv = sp.diags(r_inv)
    features = r_mat_inv.dot(features)
    return sparse_to_tuple(features)


#' in case FLAGS are defined twice
def del_all_flags(FLAGS):
    flags_dict = FLAGS._flags()
    keys_list = [keys for keys in flags_dict]
    for keys in keys_list:
        FLAGS.__delattr__(keys)

def masked_softmax_cross_entropy(preds, labels, mask):
    """Softmax cross-entropy loss with masking."""
    loss = tf.nn.softmax_cross_entropy_with_logits(logits=preds, labels=labels)
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    loss *= mask
    return tf.reduce_mean(loss)


def masked_accuracy(preds, labels, mask):
    """Accuracy with masking."""
    correct_prediction = tf.equal(tf.argmax(preds, 1), tf.argmax(labels, 1))
    accuracy_all = tf.cast(correct_prediction, tf.float32)
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    accuracy_all *= mask
    return tf.reduce_mean(accuracy_all)


def uniform(shape, scale=0.05, name=None):
    """Uniform init."""
    initial = tf.random_uniform(shape,
                                minval=-scale,
                                maxval=scale,
                                dtype=tf.float32)
    return tf.Variable(initial, name=name)


def glorot(shape, name=None):
    #Glorot深度学习中的参数初始化方法
    """Glorot & Bengio (AISTATS 2010) init."""
    init_range = np.sqrt(6.0 / (shape[0] + shape[1]))
    initial = tf.random_uniform(shape,
                                minval=-init_range,
                                maxval=init_range,
                                dtype=tf.float32)
    return tf.Variable(initial, name=name)


def zeros(shape, name=None):
    """All zeros."""
    initial = tf.zeros(shape, dtype=tf.float32)
    return tf.Variable(initial, name=name)


def ones(shape, name=None):
    """All ones."""
    initial = tf.ones(shape, dtype=tf.float32)
    return tf.Variable(initial, name=name)


def sparse_to_tuple(sparse_mx):
    """Convert sparse matrix to tuple representation."""
    def to_tuple(mx):
        if not sp.isspmatrix_coo(mx):
            mx = mx.tocoo()
        coords = np.vstack((mx.row, mx.col)).transpose()
        values = mx.data
        shape = mx.shape
        return coords, values, shape

    if isinstance(sparse_mx, list):
        for i in range(len(sparse_mx)):
            sparse_mx[i] = to_tuple(sparse_mx[i])
    else:
        sparse_mx = to_tuple(sparse_mx)

    return sparse_mx

def preprocess_adj(adj):
    """Preprocessing of adjacency matrix for scGCN model and conversion to tuple representation."""
    adj_normalized = normalize_adj(adj + sp.eye(adj.shape[0]))
    return sparse_to_tuple(adj_normalized)

def normalize_adj(adj):
    """Symmetrically normalize adjacency matrix."""
    adj = sp.coo_matrix(adj)
    rowsum = np.array(adj.sum(1))
    d_inv_sqrt = np.power(rowsum, -0.5).flatten()
    d_inv_sqrt[np.isinf(d_inv_sqrt)] = 0.
    d_mat_inv_sqrt = sp.diags(d_inv_sqrt)
    return adj.dot(d_mat_inv_sqrt).transpose().dot(d_mat_inv_sqrt).tocoo()

def construct_feed_dict(features, support, labels, labels_mask, placeholders):
    """Construct feed dictionary."""
    feed_dict = dict()
    feed_dict.update({placeholders['labels']: labels})
    feed_dict.update({placeholders['labels_mask']: labels_mask})
    feed_dict.update({placeholders['features']: features})
    feed_dict.update(
        {placeholders['support'][i]: support[i]
         for i in range(len(support))})
    feed_dict.update({placeholders['num_features_nonzero']: features[1].shape})
    return feed_dict
