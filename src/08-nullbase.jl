using Graphs,DataFrames

#  写一个函数用以生成0阶网络零模型
"""
`NullZero(net)`
Generate 0-order null model from an original network.
# Argument
* `net`:A network based on Graphs.jl.
# Return
* `null0net`:A 0-order null model network.
# Example
using Species2Networks,Graphs;
net=erdos_renyi(10, 30);
print(ne(net)==ne(NullZero(net)));
print(degree(net)==degree(NullZero(net)));
print(collect(edges(net))==collect(edges(NullZero(net))))
"""
function NullZero(net)
#  指定顶点数目和边数，生成0阶网络零模型
  null0net=Graphs.SimpleGraph(Graphs.nv(net),Graphs.ne(net))
  return(null0net)
end

#  写一个函数用以生成1阶网络零模型
"""
`NullOne(net)`
Generate 1-order null model from an original network.
# Argument
* `net`:A network based on Graphs.jl.
# Return
* `null1net`:A 1-order null model network.
# Example
using Species2Networks,Graphs;
net=erdos_renyi(10, 30);
print(ne(net)==ne(NullOne(net)));
print(degree(net)==degree(NullOne(net)));
print(collect(edges(net))==collect(edges(NullOne(net))))
"""
function NullOne(net)
#  指定节点数和度分布，生成1阶网络零模型
  null1net=Graphs.random_configuration_model(Graphs.nv(net),Graphs.degree(net))
  return(null1net)
end


"""
`NullInf(rawnet,labeln::Int,nullmode::String)`
Generate a corresponding null model network according to an original network, and obtain some properties of the null model network.
# Argument
* `rawnet`:An original network.
* `labeln`:Set the times of running label propagation algorithm.
* `nullmode`:You can choose "NullZero" or "NullOne".
# Return
* `netinf`:A dataframe includes some network properties with a null model network.
# Example
using Species2Networks,Graphs,DataFrames;
net=erdos_renyi(10, 30);
NullInf(net,100,"NullOne")
"""
function NullInf(rawnet,labeln::Int,nullmode::String)
#  对原网络进行生成对应的零模型
  if nullmode=="NullZero"
    net=NullZero(rawnet)
  elseif nullmode=="NullOne"
    net=NullOne(rawnet)
  end
#  确定需要的统计对象
  rownames=["mean betweenness centrality",
            "mean degree centrality",
            "Global clustering coefficient",
            "Average local clustering coefficient",
            "Average path distance",
            "Graph density",
            "Modularity by label propagation algorithm"]
  contentdf=DataFrame(NetworkProperties=rownames)
#  总节点数目
  tv=Graphs.nv(net)
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
  netinf=[mbc,mdc,gcc,alcc,apd,gd,md]
#  把结果数据框化
  netinfdf=DataFrame(Value=netinf)
#  合并数据框
  netinf=hcat(contentdf,netinfdf)
#  对外界返回数据
  return(netinf)
end

"""
`NullInfValue(rawnet,labeln::Int,nullmode::String)`
Generate a corresponding null model network according to an original network, and obtain some properties of the null model network.
# Argument
* `rawnet`:An original network.
* `labeln`:Set the times of running label propagation algorithm.
* `nullmode`:You can choose "NullZero" or "NullOne".
# Return
* `netinfvalue`:A dataframe includes some network properties with a null model network,but there is only value information.
# Example
using Species2Networks,Graphs,DataFrames;
net=erdos_renyi(10, 30);
NullInfValue(net,100,"NullOne")
"""
function NullInfValue(rawnet,labeln::Int,nullmode::String)
#  对原网络进行生成对应的零模型
  if nullmode=="NullZero"
    net=NullZero(rawnet)
  elseif nullmode=="NullOne"
    net=NullOne(rawnet)
  end
#  总节点数目
  tv=Graphs.nv(net)
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
  netinf=[mbc,mdc,gcc,alcc,apd,gd,md]
#  把结果数据框化
  netinfvalue=DataFrame(Value=netinf)
#  对外界返回数据
  return(netinfvalue)
end

"""
`NullTimes(rawnet,labeln::Int,nullmode::String,timen::Int)`
Generate some corresponding null model networks according to an original network, and obtain some properties of the null model networks.
# Argument
* `rawnet`:An original network.
* `labeln`:Set the times of running label propagation algorithm.
* `nullmode`:You can choose "NullZero" or "NullOne".
* `timen`:Set the times of running null model algorithm.
# Return
* `nulltimesdf`:A dataframe includes some network properties with some null model networks.
# Example
using SpeciesToNetworks,Graphs,DataFrames;
net=erdos_renyi(10, 30);
NullTimes(net,100,"NullOne",99)
"""
function NullTimes(rawnet,labeln::Int,nullmode::String,timen::Int)
#  获取第1次的网络零模型数据
  nulltimesdf=NullInf(rawnet,labeln,nullmode)
#  重新命名列名
  nulltimesdf=rename(nulltimesdf,["NetworkProperties","1"])
  for i in 2:timen
#  不断合并新的零模型数据，用string转化Int数据
    nulltimesdf=hcat(nulltimesdf,rename(NullInfValue(rawnet,labeln,nullmode),[string(i)]))
  end
#  对外界返回数据
  return(nulltimesdf)
end


