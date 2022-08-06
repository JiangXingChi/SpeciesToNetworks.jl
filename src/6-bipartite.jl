#  加载包
using MultipleTesting,HypothesisTests,StatsBase,Statistics,DataFrames,CSV

#  创建1个根据2数据框生成二分网络的相关性系数与p值的函数
"""
`BipartiteCP(dataframe1::DataFrame,dataframe2::DataFrame,method::String,adjustment::String)`
One correlation coefficient dataframe and one correlation p value dataframe are generated according to the two species dataframes,we can use these two dataframes to generate a binary network.
# Argument
* `dataframe1`:A dataframe containing species abundance information.
* `dataframe2`:A dataframe containing species abundance information.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
* `adjustment`:Select a method to adjust p value, you can use "Bonferroni","BenjaminiHochberg","BenjaminiYekutieli","BenjaminiLiu","Hochberg","Holm","Hommel","Sidak","ForwardStop","BarberCandes","raw".
# Return
* `linkcor`:A dataframe that stores the correlation coefficient.
* `linkp`:A dataframe that stores the correlation p value.
# Example
using SpeciesToNetworks,DataFrames;
a=DataFrame(a1=[1,1,0,0,0],a2=[3,3,2,2,2],a3=[1,1,2,2,2],a4=[1,2,3,4,5]);
b=DataFrame(b1=[2,2,0,0,0],b2=[3,4,2,2,1]);
linkcor,linkp=BipartiteCP(a,b,"spearman","BenjaminiHochberg");
a1b2=SpeciesCor(a[:,1],b[:,2],"spearman");
a2b1=SpeciesCor(a[:,2],b[:,1],"spearman");
a3b2=SpeciesCor(a[:,3],b[:,2],"spearman");
a4b1=SpeciesCor(a[:,4],b[:,1],"spearman");
print(linkcor);
print(a1b2);
print(a2b1);
print(a3b2);
print(a4b1)
"""
function BipartiteCP(dataframe1::DataFrame,dataframe2::DataFrame,method::String,adjustment::String)
#  获取总体物种名称
  spnames=vcat(names(dataframe1),names(dataframe2))
#  获取数据框1和数据框2的尺寸,a1与a2应当相等且为样本数量，b1与b2为两个类群的维度数目
  a1,b1=size(dataframe1)
  a2,b2=size(dataframe2)
#  计算实际需要的p值个数，用Int将其转化为整数
  cpnumber=Int(b1*b2)
#  创建空白数组array，用来储存结果
  corvalue=fill(0.0,(cpnumber))
  pvalue=fill(0.0,(cpnumber))
#  通过循环进行计算cor值和原始p值
  k=1
  for i in 1:b2
    for j in 1:b1
      corvalue[k]=SpeciesCor(dataframe2[:,i],dataframe1[:,j],method)
      pvalue[k]=SpeciesPvalue(dataframe2[:,i],dataframe1[:,j],method)
      k=k+1
    end
  end
#  对原始p值进行矫正
  adjpvalue=PvalueAdjustment(pvalue,adjustment)
#  用reshape重新包装矩阵
  rightupr=reshape(corvalue,(b1,b2))
  rightupp=reshape(adjpvalue,(b1,b2))
  leftdownr=rightupr'
  leftdownp=rightupp'
#  准备相关系数与p值矩阵的左上角与右下角
  leftupr=fill(NaN,(b1,b1))
  leftupp=fill(NaN,(b1,b1))
  rightdownr=fill(NaN,(b2,b2))
  rightdownp=fill(NaN,(b2,b2))
#  进行矩阵拼接
  linkcormatrix=vcat(hcat(leftupr,rightupr),hcat(leftdownr,rightdownr))
  linkpmatrix=vcat(hcat(leftupp,rightupp),hcat(leftdownp,rightdownp))
#  对矩阵数据框化
  linkcor=hcat(DataFrame(CorID=spnames),DataFrame(linkcormatrix,spnames))
  linkp=hcat(DataFrame(PvalueID=spnames),DataFrame(linkpmatrix,spnames))
#  对外返回数据
  return(linkcor,linkp)
end

#  创建一个根据2个物种数据框和群类名称生成网络的函数
"""
`Group2Bipartite(dataframe1::DataFrame,dataframe2::DataFrame,groupcol::Int,groupname;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector")`
The binary network was generated using two species abundance dataframes by specifying the group column and the group name.
# Argument
* `dataframe1`:A dataframe containing species abundance information.
* `dataframe2`:A dataframe containing species abundance information.
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
using SpeciesToNetworks,DataFrames;
a=DataFrame(group=["ck","ck","ck","test","test","test"],a1=[1,1,2,2,2,1],a2=[3,3,4,2,2,2],a3=[1,1,2,5,2,2],a4=[1,2,1,3,4,5]);
b=DataFrame(group=["ck","ck","ck","test","test","test"],b1=[2,2,1,1,1,2],b2=[3,4,2,2,3,1]);
indexdf,edgedf,idnetbooldf,net=Group2Bipartite(a,b,1,"ck";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector")
"""
function Group2Bipartite(dataframe1::DataFrame,dataframe2::DataFrame,groupcol::Int,groupname;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector")
#  筛选满足某个种群名称的行
  tempdf1=dataframe1[dataframe1[:,groupcol].==groupname,:]
  tempdf2=dataframe2[dataframe2[:,groupcol].==groupname,:]
#  去除种群名称的列
  tempdf1=tempdf1[:,Not(groupcol)]
   tempdf2=tempdf2[:,Not(groupcol)]
#  用条件语句处理列
  if  colfun=="RmZeroVector"
#  去除都为0的列
    tempdf1=RmZeroVector(tempdf1)
    tempdf2=RmZeroVector(tempdf2)
  elseif colfun=="AllNoZero"
#  去除丰度数据存在0的列
    tempdf1=AllNoZero(tempdf1)
    tempdf2=AllNoZero(tempdf2)
 elseif colfun=="raw"
#  不进行处理
    tempdf1=tempdf1
    tempdf2=tempdf2
  end
#  获取相关系数和相关系数p值的数据
  linkcor,linkp= BipartiteCP(tempdf1,tempdf2,method,adjustment)
#  根据相关系数和相关系数p值判断节点之间的连接情况
  indexdf,edgedf,idnetbooldf=CP2Link(linkcor,linkp,abscorrelation,pvalue)
#  根据01连接的网络矩阵转化为网络
  net=Bool2Graph(idnetbooldf)
#  对外返回数据
  return(indexdf,edgedf,idnetbooldf,net)
end

#  创建1个根据2个物种丰度数据框获取二分网络属性的函数
"""
`Groups2Bipartite(dataframe1::DataFrame,dataframe2::DataFrame,groupcol::Int;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100,writemode="NO")`
The binary network was generated using two species abundance dataframes by specifying the group column.
# Argument
* `dataframe1`:A dataframe containing species abundance information.
* `dataframe2`:A dataframe containing species abundance information.
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
using SpeciesToNetworks,DataFrames;
a=DataFrame(group=["ck","ck","ck","test","test","test"],a1=[1,1,2,2,2,1],a2=[3,3,4,2,2,2],a3=[1,1,2,5,2,2],a4=[1,2,1,3,4,5]);
b=DataFrame(group=["ck","ck","ck","test","test","test"],b1=[2,2,1,1,1,2],b2=[3,4,2,2,3,1]);
groupnetinf=Groups2Bipartite(a,b,1;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100,writemode="YES")
"""
function Groups2Bipartite(dataframe1::DataFrame,dataframe2::DataFrame,groupcol::Int;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100,writemode="NO")
#  获取种类数据
  group=levels(CategoricalArray(dataframe1[:,groupcol]))
#  获取一个初始化数据结构，用于合并之后的数据
  i=1
#  筛选满足某个种群名称的行
  tempdf1=dataframe1[dataframe1[:,groupcol].==group[i],:]
  tempdf2=dataframe2[dataframe2[:,groupcol].==group[i],:]
#  去除种群名称的列
  tempdf1=tempdf1[:,Not(groupcol)]
  tempdf2=tempdf2[:,Not(groupcol)]
#  用条件语句处理列
  if  colfun=="RmZeroVector"
#  去除都为0的列
    tempdf1=RmZeroVector(tempdf1)
    tempdf2=RmZeroVector(tempdf2)
  elseif colfun=="AllNoZero"
#  去除丰度数据存在0的列
    tempdf1=AllNoZero(tempdf1)
    tempdf2=AllNoZero(tempdf2)
  elseif colfun=="raw"
#  不进行处理
    tempdf1=tempdf1
    tempdf2=tempdf2
  end
#  获取相关系数和相关系数p值的数据
  linkcor,linkp=BipartiteCP(tempdf1,tempdf2,method,adjustment)
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
    tempdf1=dataframe1[dataframe1[:,groupcol].==group[i],:]
    tempdf2=dataframe2[dataframe2[:,groupcol].==group[i],:]
#  去除种群名称的列
    tempdf1=tempdf1[:,Not(groupcol)]
    tempdf2=tempdf2[:,Not(groupcol)]
#  用条件语句处理列
    if  colfun=="RmZeroVector"
#  去除都为0的列
      tempdf1=RmZeroVector(tempdf1)
      tempdf2=RmZeroVector(tempdf2)
    elseif colfun=="AllNoZero"
#  去除丰度数据存在0的列
      tempdf1=AllNoZero(tempdf1)
      tempdf2=AllNoZero(tempdf2)
    elseif colfun=="raw"
#  不进行处理
      tempdf1=tempdf1
      tempdf2=tempdf2
    end
#  获取相关系数和相关系数p值的数据
    linkcor,linkp=BipartiteCP(tempdf1,tempdf2,method,adjustment)
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





