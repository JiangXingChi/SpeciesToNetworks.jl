#  加载包
using CSV,DataFrames,CategoricalArrays

#  设置生成网络结构的函数
"""
`BipartiteStructure(dataframe::DataFrame,groupcol::Int;writemode="NO")`
Generate network structure.
# Argument
* `dataframe1`:A dataframe containing species abundance information.
* `dataframe2`:A dataframe containing species abundance information.
* `groupcol`:The index of the column about groups information.
# Return
*`bipartitestructure`:A dataframe containing network structure information.
# Example
using SpeciesToNetworks,DataFrames;
a=DataFrame(group=["ck","ck","ck","test","test","test"],a1=[1,1,2,2,2,1],a2=[3,3,4,2,1,2],a3=[1,1,2,5,2,2],a4=[1,2,1,3,4,5]);
b=DataFrame(group=["ck","ck","ck","test","test","test"],b1=[2,2,1,1,1,2],b2=[3,4,2,2,3,1]);
groupnetinf=Groups2Bipartite(a,b,1;adjustment="raw",abscorrelation=0,pvalue=1,writemode="YES");
bipartitestructure=BipartiteStructure(a,b,1;writemode="YES")
"""
function BipartiteStructure(dataframe1::DataFrame,dataframe2::DataFrame,groupcol::Int;writemode="NO")
#  获取节点信息
  nodename1=names(dataframe1[:,Not(groupcol)])
  nodename2=names(dataframe2[:,Not(groupcol)])
#  获取群类A个数和群类B个数
  anumber=Int(length(nodename1))
  bnumber=Int(length(nodename2))
#  获取网络一维结构的理论最大连接数目，注意用Int函数确保生成整数
  linknumber=anumber*bnumber
#  准备启动循环
#  初始化循环使用连边信息的source思想
  i=1
  source=fill(nodename1[i],bnumber)
#  进入循环，source节点的复制次数递减
  for i in 2:anumber
    tempsource=fill(nodename1[i],bnumber)
    source=vcat(source,tempsource)
  end
#  准备激动循环
#  初始化的target使用连边信息的target思想
  j=1
  target=nodename2
#  进入循环，target节点的长度从头部逐渐删除
  for j in 2:anumber
    temptarget=nodename2
    target=vcat(target,temptarget)
  end
#  使用点操作进行广播，网络一维结构的名称确定
  mergeindex=source .* "=>" .* target
#  组建一个全网络一维结构数据框
  allnetid=DataFrame(EdgeID=1:linknumber,MergeIndex=mergeindex)
#  获取群类信息
  groupinf=levels(CategoricalArray(dataframe1[:,groupcol]))
#  获取群类个数
  groupnumber=length(groupinf)
#  准备循环处理每一个类群，进行初始化
  k=1
#  设置要读取的文件名字
  tempfilename=groupinf[k]*"/"*groupinf[k]*"_edge.csv"
#  读取临时边信息
  tempedge=DataFrame(CSV.File(tempfilename))
#  使用点操作进行广播，对读取的edge信息同样生成网络一维结构名称
  partmergeindex=tempedge.SourceID .* "=>" .* tempedge.TargetID
#  组建一个局部网络一维结构数据框
  partnet=DataFrame(MergeIndex=partmergeindex,EdgeWeight=tempedge.EdgeWeight)
#  重新命名列名
  rename!(partnet,[Symbol(i) for i in ["MergeIndex",groupinf[k]]])
#  使用leftjoin对allnetid和partnet进行合并，保留allnetid的所有内容，依据MergeIndex列，生成目标群类的网络一维结构信息
  bipartitestructure=leftjoin(allnetid,partnet,on=:MergeIndex)
#  用replace对原数据进行修改
  bipartitestructure[:,k+2]=replace(bipartitestructure[:,k+2],missing=>0)
#  进入循环处理
  for k in 2:groupnumber
#  设置要读取的文件名字
    tempfilename=groupinf[k]*"/"*groupinf[k]*"_edge.csv"
#  读取临时边信息
    tempedge=DataFrame(CSV.File(tempfilename))
#  使用点操作进行广播，对读取的edge信息同样生成网络一维结构名称
    partmergeindex=tempedge.SourceID .* "=>" .* tempedge.TargetID
#  组建一个局部网络一维结构数据框
    partnet=DataFrame(MergeIndex=partmergeindex,EdgeWeight=tempedge.EdgeWeight)
#  重新命名列名
    rename!(partnet,[Symbol(i) for i in ["MergeIndex",groupinf[k]]])
#  使用leftjoin对bipartitestructure和partnet进行合并，保留bipartitestructure的所有内容，依据MergeIndex列，生成目标群类的网络一维结构信息
    bipartitestructure=leftjoin(bipartitestructure,partnet,on=:MergeIndex)
#  用replace对原数据进行修改
    bipartitestructure[:,k+2]=replace(bipartitestructure[:,k+2],missing=>0)
  end
#  用sort!在原数据上进行排序，依据EdgeID列进行排序
  sort!(bipartitestructure,[:EdgeID],rev=(false))
#  利用*符合进行string拼接，利用mkdir创建文件夹，总表写入csv文件
  if writemode=="YES"
    wd=pwd()
    bipartitestructurecsv=wd*"/bipartitestructure.csv"
    CSV.write(bipartitestructurecsv,bipartitestructure)
  end
#  对外返还数据
  return(bipartitestructure)
end