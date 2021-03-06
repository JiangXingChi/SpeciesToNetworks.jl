#  加载包
using DataFrames,Graphs

#  传递对角线数据
#  传递p值设定值
#  传递相关性设定值
"""
`CP2Link(linkcor::DataFrame,linkp::DataFrame,abscorrelation,pvalue;diagonalnumber=1.0)`
According to the correlation coefficient and P value, it is determined whether there is a link between the two species.
# Argument
* `linkcor`:A dataframe that stores the correlation coefficient.
* `linkp`:A dataframe that stores the correlation p value.
* `abscorrelation`:Set the judgment conditions of edge connection and require the minimum absolute value of correlation coefficient.
* `pvalue`:Set the judgment conditions of edge connection and require the minimum value of p value.
* `diagonalnumber`:Set the diagonal value of the network matrix. If the diagonal value is 0.0, it means that the autocorrelation of species is not considered. If the diagonal value is 1.0, it means that the autocorrelation of species is considered.
# Return
* `indexdf`:Generate index numbers for species names.
* `edgedf`:A dataframe for storing edge information.
* `idnetmatrixdf`:A dataframe for storing the network matrix.
* `idnetbooldf`:A dataframe for storing the network matrix, whether there is a connection expressed by 0 and 1.
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05)
"""
function CP2Link(linkcor::DataFrame,linkp::DataFrame,abscorrelation,pvalue;diagonalnumber=1.0)
#  去除第一列的id信息，方便ij位置的确定
  linkcorij=linkcor[:,Not(1)]
  linkpij=linkp[:,Not(1)]
#  获取linkp的尺寸，通过a判定有效物种数
  a,b=size(linkp)
#  创建ID数据
  iddf=DataFrame(ID=linkp[:,1])
#  创建索引数据
  indexdf=DataFrame(Index=1:a,ID=linkp[:,1])
#  创建一个gephi可用的矩阵，列名复制linkpij
  netmatrixdf=DataFrame(fill(diagonalnumber,(a,a)),propertynames(linkpij))
#  创建一个布尔值的初始举证，只用0，1表达连接情况
  netbooldf=DataFrame(fill(0,(a,a)),propertynames(linkpij))
#  创建一个graphs可用的数据，设置列名
  edgedf=DataFrame(Source=Int[],Target=Int[],SourceID=String[],TargetID=String[],EdgeType=String[],EdgeWeight=Float64[])
#  设置一个内循环与外循环
  for i in 1:(a-1)
    for j in (i+1):a
#  当linkpij的ij位置的p值小于等于设定的pvalue时,当linkcorij的ij位置的cor绝对值大于等于设定的abscorrelation时，判定为两个节点有连接
      if linkpij[i,j]<=pvalue && abs(linkcorij[i,j])>=abscorrelation
#  linkcorij同样的位置写入netmatrixdf，同时补写一个对称位置
        netmatrixdf[i,j]=linkcorij[i,j]
        netmatrixdf[j,i]=linkcorij[i,j]
        netbooldf[i,j]=1
        netbooldf[j,i]=1
#  设置edgedf的1、2、3、4列
        sindextemp=i
        tindextemp=j
        sidtemp=indexdf.ID[i]
        tidtemp=indexdf.ID[j]
#  对此时ij位置对应的linkcorij做判定，设置连接的type
        if linkcorij[i,j]>0.0
          ttemp="positive"
        elseif linkcorij[i,j]<0.0
          ttemp="negative"
        elseif linkcorij[i,j]==0.0
          ttemp="neutral"
        end
#  设置第6列内容，实际上就是相关系数
        wtemp=linkcorij[i,j]
#  把得到的结果组合成元组tuple，推送到目标dataframe
        tupletemp=(sindextemp,tindextemp,sidtemp,tidtemp,ttemp,wtemp)
        push!(edgedf,tupletemp)
#  当不能同时设定的p值和相关性设定值时，判定为两个两个节点没有关系，只写入0.0到netmatrixdf中，不写edgedf
      else
        netmatrixdf[i,j]=0.0
        netmatrixdf[j,i]=0.0
      end
    end
  end
#  修饰netmatrixdf
  idnetmatrixdf=hcat(iddf,netmatrixdf)
#  修饰netbooldf
  idnetbooldf=hcat(iddf,netbooldf)
#  对外返回可用于gephi和graphs的两个数据框
  return(indexdf,edgedf,idnetmatrixdf,idnetbooldf)
end

#  创建一个根据边信息和索引信息创建网络的函数
"""
`Edge2Graph(edgedf::DataFrame,indexdf::DataFrame)`
Generate a network based on a dataframe about edge information and a dataframe about vertex information.
# Argument
* `edgedf`:A dataframe about edge information.
* `indexdf`:A dataframe about vertex information.
# Return
* `net`:A network.
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Edge2Graph(edgedf,indexdf)
"""
function Edge2Graph(edgedf::DataFrame,indexdf::DataFrame)
#  获取edgedf,netmatrixdf的尺寸，a获得的就是连边数目，c获得的就是点数
  a,b=size(edgedf)
  c,d=size(indexdf)
#  创建一个无向网，确定点数为c个
  net=Graphs.SimpleGraph(c)
#  逐个添加边
  for i in 1:a
    Graphs.add_edge!(net,edgedf[i,:Source],edgedf[i,:Target])
  end
#  对外界返回网络
  return(net)
end

#  创建一个根据bool矩阵创建网络的函数
"""
`Bool2Graph(idnetbooldf::DataFrame)`
Generate a network based on a dataframe for storing the network matrix, whether there is a connection expressed by 0 and 1.
# Argument
* `idnetbooldf`:A dataframe for storing the network matrix, whether there is a connection expressed by 0 and 1.
# Return
* `net`:A network.
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf)
"""
function Bool2Graph(idnetbooldf::DataFrame)
#  去除第一列id信息后将数据框矩阵化
  netbool=Matrix(idnetbooldf[:,Not(1)])
#  根据连接矩阵转化为网络
  net=SimpleGraph(netbool)
#  对外返回网络
  return(net)
end
