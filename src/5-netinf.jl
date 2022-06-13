using DataFrames,Graphs,Statistics

#  创建网络数据统合的网络
```
`NetInf(net,edgedf::DataFrame,labeln::Int)`
Calculate some network properties by a network.
# Argument
* `net`:A network based on Graphs.jl.
* `edgedf`:A dataframe about edge information.
* `labeln`:Set the times of running label propagation algorithm.
# Return
* `netinf`:A dataframe includes some network properties.
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf);
netinf=NetInf(net,edgedf,100)
```
function NetInf(net,edgedf::DataFrame,labeln::Int)
#  确定需要的统计对象
  rownames=["Total vertices",
            "Total edges",
            "Percentage of positive edges",
            "Percentage of negative edges",
            "Average degree",
            "mean betweenness centrality",
            "mean degree centrality",
            "Global clustering coefficient",
            "Average local clustering coefficient",
            "Average path distance",
            "Graph density",
            "Modularity by label propagation algorithm"]
  contentdf=DataFrame(NetworkProperties=rownames)
#  总节点数目
  tv=Graphs.nv(net)
#  总连接数
  te=Graphs.ne(net)
#  计算正向关系占比
  ppe=(Statistics.sum(edgedf[:,:EdgeType] .== "positive")/size(edgedf,1))*100
#  计算负向关系占比
  pne=(Statistics.sum(edgedf[:,:EdgeType] .== "negative")/size(edgedf,1))*100
#  计算平均度
  ad=Statistics.mean(Graphs.degree(net))
#  计算平均介数中心性
  mbc=Statistics.mean(Graphs.betweenness_centrality(net))
#  计算平均度中心性
  mdc=Statistics.mean(Graphs.degree_centrality(net))
#  计算全局聚类系数
  gcc=Graphs.global_clustering_coefficient(net)
#  计算平均聚类系数
  alcc=Statistics.mean(Graphs.local_clustering_coefficient(net))
#  计算平均路径长度，先一个一个的点到其他点路径进行求和，在综合基础上除以2就得到实际路径和，再除以有效路径数量
  apdsumi=fill(0,tv)
  noedgei=fill(0,tv)
  for i in 1:tv
    apdtemp=Graphs.dijkstra_shortest_paths(net,i)
    apdtempd=apdtemp.dists
    noedgetemp=fill(0,tv)
    for j in 1:tv
      if apdtempd[j]>tv
        apdtempd[j]=0
        noedgetemp[j]=1
      end
    end
    apdsumi[i]=Statistics.sum(apdtempd)
    noedgei[i]=Statistics.sum(noedgetemp)
  end
  noedge=Statistics.sum(noedgei)/2
  apd=(Statistics.sum(apdsumi)/2)/((tv*(tv-1)/2)-noedge)
#  计算图密度
  gd=Graphs.density(net)
#  计算模块度，第一个参数是网络对象，第二个参数是每个节点的社团标签组成的一维array
#  执行labeln次标签传播，计算平均模块度
#  标签游走确定初始标签
  mdi=fill(0.0,labeln)
  for i in 1:labeln
    comid=Graphs.label_propagation(net)[1]
    mdi[i]=Graphs.modularity(net,comid)
  end
  md=Statistics.mean(mdi)
#  把结果组成一维array
  netinf=[tv,te,ppe,pne,ad,mbc,mdc,gcc,alcc,apd,gd,md]
#  把结果数据框化
  netinfdf=DataFrame(Value=netinf)
#  合并数据框
  netinf=hcat(contentdf,netinfdf)
#  对外界返回数据
  return(netinf)
end

#  创建网络数据值的统合函数
```
`NetInfValue(net,edgedf::DataFrame,labeln::Int)`
Calculate some network properties by a network.
# Argument
* `net`:A network based on Graphs.jl.
* `edgedf`:A dataframe about edge information.
* `labeln`:Set the times of running label propagation algorithm.
# Return
* `netinfvalue`:A dataframe includes some network properties, but there's only one column of values.
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf);
netinfvalue=NetInfValue(net,edgedf,100)
```
function NetInfValue(net,edgedf::DataFrame,labeln::Int)
#  总节点数目
  tv=Graphs.nv(net)
#  总连接数
  te=Graphs.ne(net)
#  计算正向关系占比
  ppe=(Statistics.sum(edgedf[:,:EdgeType] .== "positive")/size(edgedf,1))*100
#  计算负向关系占比
  pne=(Statistics.sum(edgedf[:,:EdgeType] .== "negative")/size(edgedf,1))*100
#  计算平均度
  ad=Statistics.mean(Graphs.degree(net))
#  计算平均介数中心性
  mbc=Statistics.mean(Graphs.betweenness_centrality(net))
#  计算平均度中心性
  mdc=Statistics.mean(Graphs.degree_centrality(net))
#  计算全局聚类系数
  gcc=Graphs.global_clustering_coefficient(net)
#  计算平均聚类系数
  alcc=Statistics.mean(Graphs.local_clustering_coefficient(net))
#  计算平均路径长度，先一个一个的点到其他点路径进行求和，在综合基础上除以2就得到实际路径和，再除以有效路径数量
  apdsumi=fill(0,tv)
  noedgei=fill(0,tv)
  for i in 1:tv
    apdtemp=Graphs.dijkstra_shortest_paths(net,i)
    apdtempd=apdtemp.dists
    noedgetemp=fill(0,tv)
    for j in 1:tv
      if apdtempd[j]>tv
        apdtempd[j]=0
        noedgetemp[j]=1
      end
    end
    apdsumi[i]=Statistics.sum(apdtempd)
    noedgei[i]=Statistics.sum(noedgetemp)
  end
  noedge=Statistics.sum(noedgei)/2
  apd=(Statistics.sum(apdsumi)/2)/((tv*(tv-1)/2)-noedge)
#  计算图密度
  gd=Graphs.density(net)
#  计算模块度，第一个参数是网络对象，第二个参数是每个节点的社团标签组成的一维array
#  执行labeln次标签传播，计算平均模块度
#  标签游走确定初始标签
  mdi=fill(0.0,labeln)
  for i in 1:labeln
    comid=Graphs.label_propagation(net)[1]
    mdi[i]=Graphs.modularity(net,comid)
  end
  md=Statistics.mean(mdi)
#  把结果组成一维array
  netinf=[tv,te,ppe,pne,ad,mbc,mdc,gcc,alcc,apd,gd,md]
#  把结果数据框化
  netinfvalue=DataFrame(Value=netinf)
#  对外界返回数据
  return(netinfvalue)
end