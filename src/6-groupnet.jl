using DataFrames,CategoricalArrays,Graphs,Statistics

#  创建一个基于丰度数据快速提取某一处理的类群，快速转化为网络
"""
`Group2Net(dataframe::DataFrame,groupcol::Int,groupname::String;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)`
Create a network based on one group.
# Argument
* `dataframe`:A dataframe containing species abundance information and the group of each sample.
* `groupcol`:The index of the column about groups information.
* `groupname`:A group name you want to study.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
* `adjustment`:Select a method to adjust p value, you can use "Bonferroni","BenjaminiHochberg","BenjaminiYekutieli","BenjaminiLiu","Hochberg","Holm","Hommel","Sidak","ForwardStop","BarberCandes","raw".
* `abscorrelation`:Set the judgment conditions of edge connection and require the minimum absolute value of correlation coefficient.
* `pvalue`:Set the judgment conditions of edge connection and require the minimum value of p value.
# Return
* `indexdf`:Generate index numbers for species names.
* `edgedf`:A dataframe for storing edge information.
* `net`:A network based on Graphs.jl.
# Example
using SpeciesToNetworks,Graphs,DataFrames,RDatasets;
dataframe=dataset("datasets", "iris");
indexdf,edgedf,net=Group2Net(dataframe,5,"setosa";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)
"""
function Group2Net(dataframe::DataFrame,groupcol::Int,groupname::String;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)
#  筛选满足某个种群名称的行
  tempdf=dataframe[dataframe[:,groupcol].==groupname,:]
#  去除种群名称的列
  tempdf=tempdf[:,Not(groupcol)]
#  去某一处理下丰度数据都为0的物种
  tempdf=RmZero(tempdf)
#  获取相关系数和相关系数p值的数据
  linkcor,linkp=SpeciesCP(tempdf,method,adjustment)
#  根据相关系数和相关系数p值判断节点之间的连接情况
  indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,abscorrelation,pvalue)
#  根据01连接的网络矩阵转化为网络
  net=Bool2Graph(idnetbooldf)
#  对外返回数据
  return(indexdf,edgedf,net)
end

#  创建一个物种丰度数据获取网络属性的函数
"""
`Groups2Netinf(dataframe::DataFrame,groupcol::Int;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,labeln=100)`
Quickly obtain the basic network information of different groups according to the species abundance dataframe.
# Argument
* `dataframe`:A dataframe containing species abundance information and the group of each sample.
* `groupcol`:The index of the column about groups information.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
* `adjustment`:Select a method to adjust p value, you can use "Bonferroni","BenjaminiHochberg","BenjaminiYekutieli","BenjaminiLiu","Hochberg","Holm","Hommel","Sidak","ForwardStop","BarberCandes","raw".
* `abscorrelation`:Set the judgment conditions of edge connection and require the minimum absolute value of correlation coefficient.
* `pvalue`:Set the judgment conditions of edge connection and require the minimum value of p value.
# Return
* `groupnetinf`:A dataframe includes basic network properties with different groups.
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("datasets", "iris");
groupnetinf=Groups2Netinf(dataframe,5;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,labeln=100)
"""
function Groups2Netinf(dataframe::DataFrame,groupcol::Int;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,labeln=100)
#  获取种类数据
  group=levels(CategoricalArray(dataframe[:,groupcol]))
#  获取一个初始化数据结构，用于合并之后的数据
  i=1
#  筛选满足某个种群名称的行
  tempdf=dataframe[dataframe[:,groupcol].==group[i],:]
#  去除种群名称的列
  tempdf=tempdf[:,Not(groupcol)]
#  去某一处理下丰度数据都为0的物种
  tempdf=RmZero(tempdf)
#  获取相关系数和相关系数p值的数据
  linkcor,linkp=SpeciesCP(tempdf,method,adjustment)
#  根据相关系数和相关系数p值判断节点之间的连接情况
  indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,abscorrelation,pvalue)
#  根据01连接的网络矩阵转化为网络
  net=Bool2Graph(idnetbooldf)
#  根据网络获取基础网络信息
  groupnetinf=NetInf(net,edgedf,labeln)
#  重现命名列名
  rename!(groupnetinf,[Symbol(i) for i in ["NetworkProperties",group[i]]])
#  从第二个处理开始循环，一个一个处理的计算网络属性
  for i in 2:size(group,1)
#  筛选满足某个种群名称的行
    tempdf=dataframe[dataframe[:,groupcol].==group[i],:]
#  去除种群名称的列
    tempdf=tempdf[:,Not(groupcol)]
#  去某一处理下丰度数据都为0的物种
    tempdf=RmZero(tempdf)
#  获取相关系数和相关系数p值的数据
    linkcor,linkp=SpeciesCP(tempdf,method,adjustment)
#  根据相关系数和相关系数p值判断节点之间的连接情况
    indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,abscorrelation,pvalue)
#  根据01连接的网络矩阵转化为网络
    net=Bool2Graph(idnetbooldf)
#  根据网络获取基础网络信息，但只是值
    tempnetinf=NetInfValue(net,edgedf,labeln)
#  重现命名列名
    rename!(tempnetinf,[Symbol(i) for i in [group[i]]])
#  合并数据框
    groupnetinf=hcat(groupnetinf,tempnetinf)
  end
#  对外返回数据
  return(groupnetinf)
end