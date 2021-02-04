import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.cross_decomposition import PLSCanonical
import random
from sklearn.neighbors import KDTree

def svd1(mat, num_cc):
    U, s, V = np.linalg.svd(mat)
    d = s[0:int(num_cc)]
    u = U[:, 0:int(num_cc)]
    v = V[0:int(num_cc), :].transpose()
    return u, v, d

def pls(x, y, num_cc):
    random.seed(42)
    plsca = PLSCanonical(n_components=int(num_cc), algorithm='svd')
    fit = plsca.fit(x, y)
    u = fit.x_weights_
    v = fit.y_weights_
    a1 = np.matmul(np.matrix(x), np.matrix(u)).transpose()
    d = np.matmul(np.matmul(a1, np.matrix(y)), np.matrix(v))
    ds = [d[i, i] for i in range(0, 30)]
    return u, v, ds


#' @import scale function
#' @rdname scale2
#' @export scaled matrix without attributes
def scale2(x):
    y = preprocessing.scale(x)
    return y

#' @param num.cc Number of canonical vectors to calculate
#' @param seed.use Random seed to set.
#' @importFrom svd1
def runcca(data1, data2, num_cc=20):
    random.seed(42)
    object1 = scale2(data1)
    object2 = scale2(data2)
    mat3 = np.matmul(np.matrix(object1).transpose(), np.matrix(object2))
    a = svd1(mat=mat3, num_cc=int(num_cc))
    cca_data = np.concatenate((a[0], a[1]))
    ind = np.where(
        [cca_data[:, col][0] < 0 for col in range(cca_data.shape[1])])[0]
    cca_data[:, ind] = cca_data[:, ind] * (-1)
    import pandas as pd
    cca_data = pd.DataFrame(cca_data)
    cca_data.index = np.concatenate(
        (np.array(data1.columns), np.array(data2.columns)))
    cca_data.columns = ['D_' + str(i) for i in range(num_cc)]
    d = a[2]
    #' d = np.around(a[2], 3)  #.astype('int')
    return cca_data, d


def l2norm(mat):
    stat = np.sqrt(np.sum(mat**2, axis=1))
    cols = mat.columns
    mat[cols] = mat[cols].div(stat, axis=0)
    mat[np.isinf(mat)] = 0
    return mat

#' @param data_use1 pandas data frame
#' @param data_use2 pandas data frame
#' @rdname runCCA
#' @export feature loadings and embeddings
def runCCA(data_use1, data_use2, features, count_names, num_cc):
    features = checkFeature(data_use1, features)
    features = checkFeature(data_use2, features)
    data1 = data_use1.loc[features, ]
    data2 = data_use2.loc[features, ]
    cca_results = runcca(data1=data1, data2=data2, num_cc=num_cc)
    cell_embeddings = np.matrix(cca_results[0])
    combined_data = data1.merge(data2,
                                left_index=True,
                                right_index=True,
                                how='inner')
    new_data1 = combined_data.loc[count_names, ].dropna()
    # loadings=loadingDim(new.data1,cell.embeddings)
    loadings = pd.DataFrame(np.matmul(np.matrix(new_data1), cell_embeddings))
    loadings.index = new_data1.index
    return cca_results, loadings


# Check if features have zero variance
# @return Returns a vector of features that is the subset of features
# that have non-zero variance
#' @param data_use pandas data frame
def checkFeature(data_use, features):
    data1 = data_use.loc[features, ]
    feature_var = data1.var(1)
    Var_features = features[np.where(feature_var != 0)[0]]
    return Var_features


#' @param data Input data
#' @param query Data to query against data
#' @param k Number of nearest neighbors to compute
def NN(data, k, query=None):
    tree = KDTree(data)
    if query is None:
        query = data
    dist, ind = tree.query(query, k)
    return dist, ind

#' @param cell_embedding : pandas data frame
def findNN(cell_embedding, cells1, cells2, k):
    print("Finding nearest neighborhoods")
    embedding_cells1 = cell_embedding.loc[cells1, ]
    embedding_cells2 = cell_embedding.loc[cells2, ]
    nnaa = NN(embedding_cells1, k=k + 1)
    nnbb = NN(embedding_cells2, k=k + 1)
    nnab = NN(data=embedding_cells2, k=k, query=embedding_cells1)
    nnba = NN(data=embedding_cells1, k=k, query=embedding_cells2)
    return nnaa, nnab, nnba, nnbb, cells1, cells2


def findMNN(neighbors, colnames, num):
    max_nn = np.array([neighbors[1][1].shape[1], neighbors[2][1].shape[1]])
    if ((num > max_nn).any()):
        num = np.min(max_nn)
        # convert cell name to neighbor index
    cells1 = colnames
    cells2 = colnames
    print("Identifying Mutual Neighbors")
    nn_cells1 = neighbors[4]
    nn_cells2 = neighbors[5]
    cell1_index = [
        list(nn_cells1).index(i) for i in cells1 if (nn_cells1 == i).any()
    ]
    cell2_index = [
        list(nn_cells2).index(i) for i in cells2 if (nn_cells2 == i).any()
    ]
    ncell = range(neighbors[1][1].shape[0])
    ncell = np.array(ncell)[np.in1d(ncell, cell1_index)]
    # initialize a list
    mnn_cell1 = [None] * (len(ncell) * num)
    mnn_cell2 = [None] * (len(ncell) * num)
    idx = -1
    for cell in ncell:
        neighbors_ab = neighbors[1][1][cell, 0:num]
        mutual_neighbors = np.where(
            neighbors[2][1][neighbors_ab, 0:num] == cell)[0]
        for i in neighbors_ab[mutual_neighbors]:
            idx = idx + 1
            mnn_cell1[idx] = cell
            mnn_cell2[idx] = i
    mnn_cell1 = mnn_cell1[0:(idx + 1)]
    mnn_cell2 = mnn_cell2[0:(idx + 1)]
    import pandas as pd
    mnns = pd.DataFrame(np.column_stack((mnn_cell1, mnn_cell2)))
    mnns.columns = ['cell1', 'cell2']
    #print("Found", mnns.shape[0], 'MNNs')
    return mnns

#' @param dim Dimension to use
#' @param numG Number of genes to return
#' @return Returns a vector of top genes
def topGenes(Loadings, dim, numG):
    data = Loadings.iloc[:, dim]
    num = np.round(numG / 2).astype('int')
    data1 = data.sort_values(ascending=False)
    data2 = data.sort_values(ascending=True)
    posG = np.array(data1.index[0:num])
    negG = np.array(data2.index[0:num])
    topG = np.concatenate((posG, negG))
    return topG


#' Get top genes across different dimensions
#' @param DimGenes How many genes to consider per dimension
#' @param maxGenes Number of genes to return at most
def TopGenes(Loadings, dims, DimGenes, maxGenes):
    maxG = max(len(dims) * 2, maxGenes)
    gens = [None] * DimGenes
    idx = -1
    for i in range(1, DimGenes + 1):
        idx = idx + 1
        selg = []
        for j in dims:
            selg.extend(set(topGenes(Loadings, dim=j, numG=i)))
        gens[idx] = set(selg)
    lens = np.array([len(i) for i in gens])
    lens = lens[lens < maxG]
    maxPer = np.where(lens == np.max(lens))[0][0] + 1
    selg = []
    for j in dims:
        selg.extend(set(topGenes(Loadings, dim=j, numG=maxPer)))
    selgene = np.array(list(set(selg)), dtype=object)
    return (selgene)


def filterPair(pairs, neighbors, mats, features, k_filter):
    nn_cells1 = neighbors[4]
    nn_cells2 = neighbors[5]
    mat1 = mats.loc[features, nn_cells1].transpose()
    mat2 = mats.loc[features, nn_cells2].transpose()
    cn_data1 = l2norm(mat1)
    cn_data2 = l2norm(mat2)
    nn = NN(data=cn_data2.loc[nn_cells2, ],
            query=cn_data1.loc[nn_cells1, ],
            k=k_filter)
    position = [
        np.where(
            pairs.loc[:, "cell2"][x] == nn[1][pairs.loc[:, 'cell1'][x], ])[0]
        for x in range(pairs.shape[0])
    ]
    nps = np.concatenate(position, axis=0)
    fpair = pairs.iloc[nps, ]
    #print("\t Finally identified ", fpair.shape[0], " MNN pairs")
    return (fpair)

def generate_graph(count_list, norm_list, scale_list, features, combine, k_filter=200, k_neighbor=5):
    all_pairs = []
    for row in combine:
        i = row[0]
        j = row[1]
        counts1 = count_list[i]
        counts2 = count_list[j]
        norm_data1 = norm_list[i]
        norm_data2 = norm_list[j]
        scale_data1 = scale_list[i]
        scale_data2 = scale_list[j]
        rowname = counts1.index
        #' @param data_use1 pandas data frame
        #' @param data_use2 pandas data frame
        #' @export feature loadings and embeddings (pandas data frame)
        cell_embedding, loading = runCCA(data_use1=scale_data1,
                                         data_use2=scale_data2,
                                         features=features,
                                         count_names=rowname,
                                         num_cc=30)
        norm_embedding = l2norm(mat=cell_embedding[0])
        #' identify nearest neighbor
        cells1 = counts1.columns
        cells2 = counts2.columns
        neighbor = findNN(cell_embedding=norm_embedding,
                          cells1=cells1,
                          cells2=cells2,
                          k=30)
        #' identify mutual nearest neighbors
        #' @param neighbors,colnames
        #' @export mnn_pairs
        mnn_pairs = findMNN(neighbors=neighbor,
                            colnames=cell_embedding[0].index,
                            num=k_neighbor)
        select_genes = TopGenes(Loadings=loading,
                                dims=range(30),
                                DimGenes=100,
                                maxGenes=200)
        Mat = pd.concat([norm_data1, norm_data2], axis=1)
        final_pairs = filterPair(pairs=mnn_pairs,
                                 neighbors=neighbor,
                                 mats=Mat,
                                 features=select_genes,
                                 k_filter=k_filter)
        final_pairs['Dataset1'] = [i + 1] * final_pairs.shape[0]
        final_pairs['Dataset2'] = [j + 1] * final_pairs.shape[0]
        all_pairs.append(final_pairs)
    return all_pairs
