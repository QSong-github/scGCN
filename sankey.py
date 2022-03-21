import pandas as pd
from pyecharts.charts import Sankey
from pyecharts import options as opts
# pd.set_option('display.max_columns',20)
# pd.set_option('display.max.rows',20000)

def createSankey(dir_in,rel_species="original",qur_species="predicted",species="species"):
    # data_in=pd.read_csv(dir_in,header=None)
    data_in=dir_in
    data=pd.DataFrame(columns=["original","predicted","number"])
    aa=data_in["predicted"].groupby(data_in["original"]).value_counts()
    # print(aa)
    cc=[]
    dd=[]
    for i in range(len(aa.index)):
        str_a=str(aa.index[i][0])
        cc.append(str_a)
        str_b=str(aa.index[i][1])
        dd.append(str_b)
    bb=aa.values.tolist()
    data["original"]=cc
    data["predicted"]=dd
    data["number"]=bb
    # print(data_out)

    data.columns=["rel","pre","num"]
    data.dropna(axis=0,how="any",inplace=True) #生成一个没有空值的一对一同源的pd表格
    # print(data.head())
    def changeNames1(dataa):
        a=[]
        a= rel_species+"_"+dataa['rel']
        return a
    def changeNames2(dataa):
        a=[]
        a=qur_species+"_"+dataa['pre']
        return a
        # dataa["num"]=dataa["num"]
    print("not ok")
    data["rel"]=data.apply(lambda r:changeNames1(r),axis=1)
    data["pre"]=data.apply(lambda r:changeNames2(r),axis=1)
    ###################################################################
    def createNode(a,b):
        nodes=[]
        for i in range(0, len(a)):
            dic = {}
            if(len(a[i])!=0):

                dic["name"] = a[i]
                nodes.append(dic)
            else:
                pass
        for j in range(0, len(b)):
            dic={}
            if(len(b[j])!=0):

                dic["name"] = b[j]
                nodes.append(dic)
            else:
                pass
        return nodes
    sig_rel=list(set(data["rel"].tolist()))
    sig_pre=list(set(data["pre"].tolist()))
    # sig_rel=sig_data["rel"].tolist()
    # sig_pre=sig_data["pre"].tolist()
    sig_nodes=createNode(sig_rel,sig_pre)
    # print(sig_nodes)
    def createLinks(dataa):
        links = []
        for i in dataa.values:
            dic = {}
            if(len(i[0])!=0):
                dic["source"]=i[0]
                dic["target"]=i[1]
                dic["value"]=i[2]
                links.append(dic)
            else:
                pass
        return links
    sig_links=createLinks(data)
    # print(sig_links)
    ################################################################
    c=Sankey()
    c.add(
        series_name="the number of cells",
        nodes=sig_nodes,
        links=sig_links,
        pos_top="10%",  #设置图片离标题的高度
        #pos_right="5%",
        #pos_bottom="5%",
        pos_left="5%",
        is_draggable=False, #设置是否可以拖动节点
        focus_node_adjacency=True,
        linestyle_opt=opts.LineStyleOpts(opacity=0.2, curve=0.5, color="source",type_="dotted"),
        label_opts=opts.LabelOpts(position="right",),
    )
    c.set_global_opts(
        title_opts=opts.TitleOpts(title="left:"+rel_species+";"+"right:"+qur_species),
        # toolbox_opts=opts.ToolboxOpts(
        #     pos_right="20%",
        #)
        # legend_opts=opts.LegendOpts(
        #     is_show=True,
        # )
        # tooltip_opts=opts.TooltipOpts(
        #
        # ),
        # visualmap_opts=opts.VisualMapOpts(
        #     is_calculable=True,
        # )
        # datazoom_opts=opts.DataZoomOpts(
        #     is_show=True,
        # )

    )
    # c = (
    #         Sankey()
    #         .add(
    #             series_name="the number of cells",
    #             nodes=sig_nodes,
    #             links=sig_links,
    #             linestyle_opt=opts.LineStyleOpts(opacity=0.2, curve=0.5, color="source",type_="dotted"),
    #             label_opts=opts.LabelOpts(position="right",),
    #         )
    #         .set_global_opts(title_opts=opts.TitleOpts(title="left:"+rel_species+";"+"right:"+qur_species))
    #     )
    # 输出html可视化结果
    dir_out="./"+species+"_sankey.html"
    c.render(dir_out)
# createGraph("./datain/data_in.txt")