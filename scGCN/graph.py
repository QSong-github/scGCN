from utility import *
import glob
import pandas as pd
import os

def graph_construct(outputdir):
    path0 = os.path.join(os.getcwd(), outputdir)

    #' import processed data 
    files1 = glob.glob(path0 + "/count_data/*.csv")
    files1.sort()
    count_list = []
    for df in files1:
        print(df)
        count_list.append(pd.read_csv(df, index_col=0))

    files2 = glob.glob(path0 + "/norm_data/*.csv")
    files2.sort()
    norm_list = []
    for df in files2:
        print(df)
        norm_list.append(pd.read_csv(df, index_col=0))

    files3 = glob.glob(path0 + "/scale_data/*.csv")
    files3.sort()
    scale_list = []
    for df in files3:
        print(df)
        scale_list.append(pd.read_csv(df, index_col=0))

    files4 = glob.glob(path0 + "/label_data/*.csv")
    files4.sort()
    label_list = []
    for df in files4:
        print(df)
        label_list.append(pd.read_csv(df, index_col=0))

    fpath = os.path.join(path0, 'sel_features.csv')
    features = pd.read_csv(fpath, index_col=False).values.flatten()

    #' graph construction    
    import itertools

    N = len(count_list)
    if (N == 1):
        combine = pd.Series([(0, 0)])
    else:
        combin = list(itertools.product(list(range(N)), list(range(N))))
        index = [i for i, x in enumerate([i[0] < i[1] for i in combin]) if x]
        combine = pd.Series(combin)[index]

    pairss1 = generate_graph(count_list=count_list,
                             norm_list=norm_list,
                             scale_list=scale_list,
                             features=features,
                             combine=combine, k_neighbor=10)

    count_list2 = [count_list[1], count_list[1]]
    norm_list2 = [norm_list[1], norm_list[1]]
    scale_list2 = [scale_list[1], scale_list[1]]

    pairss2 = generate_graph(count_list=count_list2,
                             norm_list=norm_list2,
                             scale_list=scale_list2,
                             features=features,
                             combine=combine,k_neighbor=10)

    #'@param graph1: inter-dataset graph 

    graph1 = pairss1[0].iloc[:, 0:2].reset_index()
    graph1.to_csv('./input/inter_graph.csv')

    #'@param graph2: intra-dataset graph 
    graph2 = pairss2[0].iloc[:, 0:2].reset_index()
    graph2.to_csv('./input/intra_graph.csv')

    label1 = label_list[0]
    label1.to_csv('./input/Label1.csv', index=False)

    label2 = label_list[1]
    label2.to_csv('./input/Label2.csv', index=False)
