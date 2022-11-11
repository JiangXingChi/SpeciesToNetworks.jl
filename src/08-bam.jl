#  创建1个根据2数据框生成二分网络的相关性系数与p值的函数
"""
`BipartiteAdjCP(dataframe1::DataFrame,dataframe2::DataFrame,method::String,adjustment::String)`
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
linkcor,linkp=BipartiteAdjCP(a,b,"spearman","BenjaminiHochberg");
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
function BipartiteAdjCP(dataframe1::DataFrame,dataframe2::DataFrame,method::String,adjustment::String)
#  获取物种名称
  df1names=names(dataframe1)
  df2names=names(dataframe2)
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
#  对矩阵数据框化
  linkcor=hcat(DataFrame(CorID=df1names),DataFrame(rightupr,df2names))
  linkp=hcat(DataFrame(PvalueID=df1names),DataFrame(rightupp,df2names))
#  对外返回数据
  return(linkcor,linkp)
end

#  传递对角线数据
#  传递p值设定值
#  传递相关性设定值
"""
`BACP2Link(linkcor::DataFrame,linkp::DataFrame,abscorrelation,pvalue)`
According to the correlation coefficient and P value, it is determined whether there is a link between the two species.
# Argument
* `linkcor`:A dataframe that stores the correlation coefficient.
* `linkp`:A dataframe that stores the correlation p value.
* `abscorrelation`:Set the judgment conditions of edge connection and require the minimum absolute value of correlation coefficient.
* `pvalue`:Set the judgment conditions of edge connection and require the maximum value of p value.
# Return
* `indexdf`:Generate index numbers for species names.
* `edgedf`:A dataframe for storing edge information.
* `idnetbooldf`:A dataframe for storing the adjacency matrix, whether there is a connection expressed by 0(false) and 1(true).
# Example
using SpeciesToNetworks,DataFrames;
a=DataFrame(a1=[1,1,0,0,0],a2=[3,3,2,2,2],a3=[1,1,2,2,2],a4=[1,2,3,4,5]);
b=DataFrame(b1=[2,2,0,0,0],b2=[3,4,2,2,1]);
linkcor,linkp=BipartiteAdjCP(a,b,"spearman","BenjaminiHochberg");
print(linkcor,linkp);
indexdf,edgedf,idnetbooldf=BACP2Link(linkcor,linkp,0.6,0.05)
"""
function BACP2Link(linkcor::DataFrame,linkp::DataFrame,abscorrelation,pvalue)
#  去除第一列的id信息，方便ij位置的确定
  linkcorij=linkcor[:,Not(1)]
  linkpij=linkp[:,Not(1)]
#  获取linkpij的尺寸，通过a判定有效物种数
  a,b=size(linkpij)
#  创建ID数据
  mn1=linkp[:,1]
  mn2=names(linkpij)
  membername=vcat(mn1,mn2)
  groupclass=vcat(fill("GroupA",a),fill("GroupB",b))
  iddf=DataFrame(Index=1:(a+b),ID=membername,Group=groupclass)
#  创建一个布尔值的初始邻接矩阵，只用0，1表达连接情况
  netbooldf=DataFrame(fill(0,(a,b)),propertynames(linkpij))
#  创建一个graphs可用的数据，设置列名
  edgedf=DataFrame(SourceIndex=Int[],TargetIndex=Int[],Source=String[],Target=String[],EdgeType=String[],EdgeWeight=Float64[])
#  设置一个内循环与外循环
  for i in 1:a
    for j in 1:b
#  当linkpij的ij位置的p值小于等于设定的pvalue时,当linkcorij的ij位置的cor绝对值大于等于设定的abscorrelation时，判定为两个节点有连接
      if linkpij[i,j]<=pvalue && abs(linkcorij[i,j])>=abscorrelation
#  将1值（true）写入邻接矩阵，同时补写一个对称位置
        netbooldf[i,j]=1
#  设置edgedf的1、2、3、4列
        sindextemp=i
        tindextemp=j+a
        sidtemp=mn1[i]
        tidtemp=mn2[j]
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
      end
    end
  end
#  修饰netbooldf
  idnetbooldf=hcat(DataFrame(BipartiteAdjacencyMatrix=mn1),netbooldf)
#  对外返回可用于gephi和graphs的两个数据框
  return(iddf,edgedf,idnetbooldf)
end

#  创建一个根据2个物种数据框和群类名称生成网络的函数
"""
`Group2BAM(dataframe1::DataFrame,dataframe2::DataFrame,groupcol::Int,groupname;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector")`
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
# Example
using SpeciesToNetworks,DataFrames;
a=DataFrame(group=["ck","ck","ck","test","test","test"],a1=[1,1,2,2,2,1],a2=[3,3,4,2,2,2],a3=[1,1,2,5,2,2],a4=[1,2,1,3,4,5]);
b=DataFrame(group=["ck","ck","ck","test","test","test"],b1=[2,2,1,1,1,2],b2=[3,4,2,2,3,1]);
indexdf,edgedf,idnetbooldf=Group2BAM(a,b,1,"ck";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector")
"""
function Group2BAM(dataframe1::DataFrame,dataframe2::DataFrame,groupcol::Int,groupname;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector")
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
  linkcor,linkp=BipartiteAdjCP(tempdf1,tempdf2,method,adjustment)
#  根据相关系数和相关系数p值判断节点之间的连接情况
  indexdf,edgedf,idnetbooldf=BACP2Link(linkcor,linkp,abscorrelation,pvalue)
#  对外返回数据
  return(indexdf,edgedf,idnetbooldf)
end

#  创建1个根据2个物种丰度数据框获取二分网络属性的函数
"""
`Groups2BAM(dataframe1::DataFrame,dataframe2::DataFrame,groupcol::Int;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100,writemode="NO")`
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
* `writemode`:Whether to write out the point, edge and adjacency data of the group class as CSV files, with "YES" and "NO" modes.
# Return
* `indexdf`:Generate index numbers for species names.
* `edgedf`:A dataframe for storing edge information.
* `idnetbooldf`:A dataframe for storing the adjacency matrix, whether there is a connection expressed by 0(false) and 1(true).
# Example
using SpeciesToNetworks,DataFrames;
a=DataFrame(group=["ck","ck","ck","test","test","test"],a1=[1,1,2,2,2,1],a2=[3,3,4,2,2,2],a3=[1,1,2,5,2,2],a4=[1,2,1,3,4,5]);
b=DataFrame(group=["ck","ck","ck","test","test","test"],b1=[2,2,1,1,1,2],b2=[3,4,2,2,3,1]);
groupnetinf=Groups2BAM(a,b,1;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",writemode="YES")
"""
function Groups2BAM(dataframe1::DataFrame,dataframe2::DataFrame,groupcol::Int;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",writemode="NO")
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
  linkcor,linkp=BipartiteAdjCP(tempdf1,tempdf2,method,adjustment)
#  根据相关系数和相关系数p值判断节点之间的连接情况
  indexdf,edgedf,idnetbooldf=BACP2Link(linkcor,linkp,abscorrelation,pvalue)
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
  linkcor,linkp=BipartiteAdjCP(tempdf1,tempdf2,method,adjustment)
#  根据相关系数和相关系数p值判断节点之间的连接情况
  indexdf,edgedf,idnetbooldf=BACP2Link(linkcor,linkp,abscorrelation,pvalue)
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
#  对外返回数据
  return(indexdf,edgedf,idnetbooldf)
end