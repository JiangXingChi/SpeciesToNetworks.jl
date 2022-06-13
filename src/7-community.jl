using Graphs,DataFrames,Statistics,CategoricalArrays

#  设置一个逐步合并社团的函数
```
`ReduceLabel(comid::Vector,net)`
Select two communities to merge, the evaluation standard is modularity Q value.
# Argument
* `comid`:Original vertex label.
* `net`:A network based on Graphs.jl.
# Return
* `newcomid`:New vertex label.
* `Q`:New modularity Q value.
* `qvalue`:Old modularity Q value.
* `newcomnumber`:New community numbers.
* `comnumber`:Old community numbers.
# Example
using SpeciesToNetworks,DataFrames,RDatasets,Graphs;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf);
comid=Graphs.label_propagation(net)[1];
newcomid,Q,qvalue,newcomnumber,comnumber=ReduceLabel(comid,net)
```
function ReduceLabel(comid::Vector,net)
#  获取初始社团数目
  comnumber=findmax(comid)[1]
#  计算此时的模块Q值
  qvalue=Graphs.modularity(net,comid)
#  创建一个社团数平方的矩阵，用以储存两两合并的新结果，初始值为初始标签传播的Q值
  newqmatrix=fill(0.0,(comnumber,comnumber))
#  用ij双循环处理两两社团合并
  for i in 1:comnumber-1
    for j in i+1:comnumber
#  把j替换成i，组成一个新的临时标签向量
      tempcomid=replace(comid,j=>i)
#  基于新的临时标签向量计算新的临时模块度Q
      tempqvalue=Graphs.modularity(net,tempcomid)
#  把新的临时模块度Q填进newqmatrix
      newqmatrix[i,j]=tempqvalue
      newqmatrix[j,i]=tempqvalue
    end
  end
#  确定两两合并后的最大Q值，用[1]锁定最大值
    Q=findmax(newqmatrix)[1]
#  当某两两合并的新Q值大于等于一开始的Q值时，确认启动合并
  if Q >= qvalue
#  用findmax找到最大Q值，用[2]锁定笛卡尔坐标，用[1][2]分别提取笛卡尔的一维二维坐标值
    a=findmax(newqmatrix)[2][1]
    b=findmax(newqmatrix)[2][2]
#  确定新的标签序列，把大的社团标签转换为小的社团标签，这里的replace不加!，不直接影响
    newcomid=replace(comid,max(a,b)=>min(a,b))
#  把大的空出的社团标签后的标签一次减1
    for k in (max(a,b)+1):comnumber
#  这里的replace加!，直接影响
      replace!(newcomid,k=>(k-1))
    end
#  计算此时的社团数目
    newcomnumber=findmax(newcomid)[1]
#  当合并后的Q值不如原来Q值时，新标签向量复制原标签向量就可以了
  else
    newcomid=copy(comid)
#  确定此时的社团数目
    newcomnumber=copy(comnumber)
#  确定此时的Q值
    Q=copy(qvalue)
  end
#  返回新标签信息向量，新的模块度Q值，旧的模块度Q值，新的社团数量，旧的社团数量
  return(newcomid,Q,qvalue,newcomnumber,comnumber)
end

#  设置一个社团合并到极限的函数
```
`MinLabel(comid::Vector,net)`
Select maximum communities to merge, the evaluation standard is modularity Q value.
# Argument
* `comid`:Original vertex label.
* `net`:A network based on Graphs.jl.
# Return
* `newcomid`:New vertex label.
* `Q`:New modularity Q value.
* `beginqvalue`:Initial modularity Q value.
* `newcomnumber`:New community numbers.
* `begincomnumber`:Initial community numbers.
# Example
using SpeciesToNetworks,DataFrames,RDatasets,Graphs;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf);
comid=Graphs.label_propagation(net)[1];
newcomid,Q,beginqvalue,newcomnumber,begincomnumber=MinLabel(comid,net)
```
function MinLabel(comid::Vector,net)
#  计算初始的模块度Q值
  beginqvalue=Graphs.modularity(net,comid)
#  计算初始的社团数目
  begincomnumber=findmax(comid)[1]
#  执行第一次社团合并，获得相应的值
  newcomid,Q,qvalue,newcomnumber,comnumber=ReduceLabel(comid,net)
#  当新的社团数量小于旧的社团数量时，持续执行合并，否则中断
  while newcomnumber < comnumber
#  执行新的合并之后，更新一系列值，这里主要更新新的社团数量和旧的社团数量
    newcomid,Q,qvalue,newcomnumber,comnumber=ReduceLabel(newcomid,net)
  end
#  合并到极限后，返回新的标签向量，新的模块度Q值，初始的模块度Q值，新的社团数量，初始的社团数量
  return(newcomid,Q,beginqvalue,newcomnumber,begincomnumber)
end

#  设置一个获取顶点标签数据的函数
```
`LabelDf(dataframe::DataFrame,groupcol::Int,groupname::String;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)`
Create a label information based on one group.
# Argument
* `dataframe`:A dataframe containing species abundance information and the group of each sample.
* `groupcol`:The index of the column about groups information.
* `groupname`:A group name you want to study.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
* `adjustment`:Select a method to adjust p value, you can use "Bonferroni","BenjaminiHochberg","BenjaminiYekutieli","BenjaminiLiu","Hochberg","Holm","Hommel","Sidak","ForwardStop","BarberCandes","raw".
* `abscorrelation`:Set the judgment conditions of edge connection and require the minimum absolute value of correlation coefficient.
* `pvalue`:Set the judgment conditions of edge connection and require the minimum value of p value.
# Return
* `labeldf`:A dataframe containing vertex index, species information and two types of label information.
* `edgedf`:A dataframe for storing edge information.
* `Q`:New modularity Q value.
* `beginqvalue`:Initial modularity Q value.
* `newcomnumber`:New community numbers.
* `begincomnumber`:Initial community numbers.
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("datasets", "iris");
labeldf,edgedf,Q,beginqvalue,newcomnumber,begincomnumber=LabelDf(dataframe,5,"setosa";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)
```
function LabelDf(dataframe::DataFrame,groupcol::Int,groupname::String;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)
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
#  利用标签传播法获取初始标签序列
  comid=Graphs.label_propagation(net)[1]
#  把初始标签序列整合成数据框
  comiddf=DataFrame(LabelPropagation=comid)
#  用MinLable将社团两两合并到极限
  newcomid,Q,beginqvalue,newcomnumber,begincomnumber=MinLabel(comid,net)
#  将合并到极限的标签序列整合成成数据框
  newcomiddf=DataFrame(LabelMaxQ=newcomid)
#  合并数据框
  labeldf=hcat(indexdf,comiddf,newcomiddf)
#  对外部返回整合的数据框，新的模块度Q值，初始的模块度Q值，新的社团数量，初始的社团数量
  return(labeldf,edgedf,Q,beginqvalue,newcomnumber,begincomnumber)
end
