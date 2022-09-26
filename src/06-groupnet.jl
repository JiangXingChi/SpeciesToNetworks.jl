using DataFrames,CategoricalArrays,Graphs,Statistics,CSV

#  创建一个基于丰度数据快速提取某一处理的类群，快速转化为网络
"""
`Group2Net(dataframe::DataFrame,groupcol::Int,groupname;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector")`
Create a network based on one group.
# Argument
* `dataframe`:A dataframe containing species abundance information and sample groups information.
* `groupcol`:The index of the column about groups information.
* `groupname`:A group name you want to study.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
* `adjustment`:Select a method to adjust p value, you can use "Bonferroni","BenjaminiHochberg","BenjaminiYekutieli","BenjaminiLiu","Hochberg","Holm","Hommel","Sidak","ForwardStop","BarberCandes","raw".
* `abscorrelation`:Set the judgment conditions of edge connection and require the minimum absolute value of correlation coefficient.
* `pvalue`:Set the judgment conditions of edge connection and require the maximum value of p value.
* `colfun`:Set the function to process the column,you can use "RmZeroVector","AllNoZero","raw".
# Return
* `indexdf`:Generate index numbers for species names.
* `edgedf`:A dataframe for storing edge information.
* `idnetbooldf`:A dataframe for storing the adjacency matrix, whether there is a connection expressed by 0(false) and 1(true).
* `net`:A network based on Graphs.jl.
# Example
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(Groups=["a","a","a","a","b","b","b","b"],species1=[2,2,1,1,0,5,7,2],species2=[0,0,0,0,3,2,2,2],species3=[1,1,2,2,6,8,2,2],species4=[0,2,2,4,9,3,4,5]);
indexdf1,edgedf1,idnetbooldf1,net1=Group2Net(dataframe,1,"a";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector");
indexdf2,edgedf2,idnetbooldf2,net2=Group2Net(dataframe,1,"a";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="AllNoZero");
indexdf3,edgedf3,idnetbooldf3,net3=Group2Net(dataframe,1,"a";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="raw");
print(idnetbooldf1,idnetbooldf2,idnetbooldf3)
"""
function Group2Net(dataframe::DataFrame,groupcol::Int,groupname;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector")
#  筛选满足某个种群名称的行
  tempdf=dataframe[dataframe[:,groupcol].==groupname,:]
#  去除种群名称的列
  tempdf=tempdf[:,Not(groupcol)]
#  用条件语句处理列
  if  colfun=="RmZeroVector"
#  去除都为0的列
    tempdf=RmZeroVector(tempdf)
  elseif colfun=="AllNoZero"
#  去除丰度数据存在0的列
    tempdf=AllNoZero(tempdf)
  elseif colfun=="raw"
#  不进行处理
    tempdf=tempdf
  end
#  获取相关系数和相关系数p值的数据
  linkcor,linkp=SpeciesCP(tempdf,method,adjustment)
#  根据相关系数和相关系数p值判断节点之间的连接情况
  indexdf,edgedf,idnetbooldf=CP2Link(linkcor,linkp,abscorrelation,pvalue)
#  根据01连接的网络矩阵转化为网络
  net=Bool2Graph(idnetbooldf)
#  对外返回数据
  return(indexdf,edgedf,idnetbooldf,net)
end

#  创建一个物种丰度数据获取网络属性的函数
"""
`Groups2Net(dataframe::DataFrame,groupcol::Int;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100, writemode="NO")`
Quickly obtain the basic network information of different groups according to the species abundance dataframe.
# Argument
* `dataframe`:A dataframe containing species abundance information and the group of each sample.
* `groupcol`:The index of the column about groups information.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
* `adjustment`:Select a method to adjust p value, you can use "Bonferroni","BenjaminiHochberg","BenjaminiYekutieli","BenjaminiLiu","Hochberg","Holm","Hommel","Sidak","ForwardStop","BarberCandes","raw".
* `abscorrelation`:Set the judgment conditions of edge connection and require the minimum absolute value of correlation coefficient.
* `pvalue`:Set the judgment conditions of edge connection and require the maximum value of p value.
* `colfun`:Set the function to process the column,you can use "RmZeroVector","AllNoZero","raw".
* `labeln`:Set the times of running label propagation algorithm.
* `writemode`:Whether to write out the point, edge and adjacency data of the group class as CSV files, with "YES" and "NO" modes.
# Return
* `groupnetinf`:A dataframe includes basic network properties with different groups.
# Example
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(Groups=["a","a","a","a","b","b","b","b"],species1=[2,2,1,1,0,5,7,2],species2=[0,0,0,0,3,2,2,2],species3=[1,1,2,2,6,8,2,2],species4=[0,2,2,4,9,3,4,5]);
groupnetinf=Groups2Net(dataframe,1;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100,writemode="YES")
"""
function Groups2Net(dataframe::DataFrame,groupcol::Int;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100,writemode="NO")
#  获取种类数据
  group=levels(CategoricalArray(dataframe[:,groupcol]))
#  获取一个初始化数据结构，用于合并之后的数据
  i=1
#  筛选满足某个种群名称的行
  tempdf=dataframe[dataframe[:,groupcol].==group[i],:]
#  去除种群名称的列
  tempdf=tempdf[:,Not(groupcol)]
#  用条件语句处理列
  if  colfun=="RmZeroVector"
#  去除都为0的列
    tempdf=RmZeroVector(tempdf)
  elseif colfun=="AllNoZero"
#  去除丰度数据存在0的列
    tempdf=AllNoZero(tempdf)
  elseif colfun=="raw"
#  不进行处理
    tempdf=tempdf
  end
#  获取相关系数和相关系数p值的数据
  linkcor,linkp=SpeciesCP(tempdf,method,adjustment)
#  根据相关系数和相关系数p值判断节点之间的连接情况
  indexdf,edgedf,idnetbooldf=CP2Link(linkcor,linkp,abscorrelation,pvalue)
#  根据01连接的网络矩阵转化为网络
  net=Bool2Graph(idnetbooldf)
#  根据网络获取基础网络信息
  groupnetinf=NetInf(net,edgedf,labeln)
#  重现命名列名
  rename!(groupnetinf,[Symbol(i) for i in ["NetworkProperties",group[i]]])
#  利用*符合进行string拼接，利用mkdir创建文件夹，按群类写入csv文件
  if writemode=="YES"
    wd=pwd()
    mkdir(group[i])
     indexdfcsv=wd*"/"*group[i]*"/"*group[i]*"_node"*".csv"
     edgedfcsv=wd*"/"*group[i]*"/"*group[i]*"_edge"*".csv"
     idnetbooldfcsv=wd*"/"*group[i]*"/"*group[i]*"_adjacency"*".csv"
     CSV.write( indexdfcsv,indexdf)
     CSV.write( edgedfcsv,edgedf)
     CSV.write( idnetbooldfcsv,idnetbooldf)
  end
#  从第二个处理开始循环，一个一个处理的计算网络属性
  for i in 2:size(group,1)
#  筛选满足某个种群名称的行
    tempdf=dataframe[dataframe[:,groupcol].==group[i],:]
#  去除种群名称的列
    tempdf=tempdf[:,Not(groupcol)]
#  用条件语句处理列
    if  colfun=="RmZeroVector"
#  去除都为0的列
      tempdf=RmZeroVector(tempdf)
    elseif colfun=="AllNoZero"
#  去除丰度数据存在0的列
      tempdf=AllNoZero(tempdf)
    elseif colfun=="raw"
#  不进行处理
      tempdf=tempdf
    end
#  获取相关系数和相关系数p值的数据
    linkcor,linkp=SpeciesCP(tempdf,method,adjustment)
#  根据相关系数和相关系数p值判断节点之间的连接情况
    indexdf,edgedf,idnetbooldf=CP2Link(linkcor,linkp,abscorrelation,pvalue)
#  根据01连接的网络矩阵转化为网络
    net=Bool2Graph(idnetbooldf)
#  根据网络获取基础网络信息，但只是值
    tempnetinf=NetInfValue(net,edgedf,labeln)
#  重现命名列名
    rename!(tempnetinf,[Symbol(i) for i in [group[i]]])
#  合并数据框
    groupnetinf=hcat(groupnetinf,tempnetinf)
#  利用*符合进行string拼接，利用mkdir创建文件夹，按群类写入csv文件
    if writemode=="YES"
      wd=pwd()
      mkdir(group[i])
      indexdfcsv=wd*"/"*group[i]*"/"*group[i]*"_node"*".csv"
      edgedfcsv=wd*"/"*group[i]*"/"*group[i]*"_edge"*".csv"
      idnetbooldfcsv=wd*"/"*group[i]*"/"*group[i]*"_adjacency"*".csv"
      CSV.write( indexdfcsv,indexdf)
      CSV.write( edgedfcsv,edgedf)
      CSV.write( idnetbooldfcsv,idnetbooldf)
  end
  end
#  利用*符合进行string拼接，利用mkdir创建文件夹，总表写入csv文件
  if writemode=="YES"
     wd=pwd()
     groupnetinfcsv=wd*"/GroupsInformation.csv"
     CSV.write( groupnetinfcsv,groupnetinf)
  end
#  对外返回数据
  return(groupnetinf)
end
