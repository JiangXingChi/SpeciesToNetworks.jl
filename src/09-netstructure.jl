#  加载包
using CSV,DataFrames,CategoricalArrays

#  设置生成网络结构的函数
"""
`NetStructure(dataframe::DataFrame,groupcol::Int;writemode="NO")`
Generate network structure.
# Argument
* `dataframe`:A dataframe containing species abundance information and sample groups information.
* `groupcol`:The index of the column about groups information.
# Return
*`netstructure`:A dataframe containing network structure information.
# Example
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(Groups=["a","a","a","a","b","b","b","b"],species1=[2,2,1,1,0,5,7,2],species2=[0,0,0,0,3,2,2,2],species3=[1,1,2,2,6,8,2,2],species4=[0,2,2,4,9,3,4,5]);
groupnetinf=Groups2Net(dataframe,1;abscorrelation=0,pvalue=1,writemode="YES");
netstructure=NetStructure(dataframe,1;writemode="YES")
"""
function NetStructure(dataframe::DataFrame,groupcol::Int;writemode="NO")
#  获取节点信息
  nodename=names(dataframe[:,Not(groupcol)])
#  获取节点数目
  nodenumber=Int(length(nodename))
#  获取网络一维结构的理论最大连接数目，注意用Int函数确保生成整数
  linknumber=Int(nodenumber*(nodenumber-1)/2)
#  准备启动循环
#  初始化循环使用连边信息的source思想
  i=1
  source=fill(nodename[i],(nodenumber-i))
#  进入循环，source节点的复制次数递减
  for i in 2:(nodenumber-1)
    tempsource=fill(nodename[i],(nodenumber-i))
    source=vcat(source,tempsource)
  end
#  准备激动循环
#  初始化的target使用连边信息的target思想
  j=1
  target=nodename[(j+1):nodenumber]
#  进入循环，target节点的长度从头部逐渐删除
  for j in 2:(nodenumber-1)
    temptarget=nodename[(j+1):nodenumber]
    target=vcat(target,temptarget)
  end
#  使用点操作进行广播，网络一维结构的名称确定
  mergeindex=source .* "=>" .* target
#  组建一个全网络一维结构数据框
  allnetid=DataFrame(EdgeID=1:linknumber,MergeIndex=mergeindex)
#  获取群类信息
  groupinf=levels(CategoricalArray(dataframe[:,groupcol]))
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
  netstructure=leftjoin(allnetid,partnet,on=:MergeIndex)
#  用replace对原数据进行修改
  netstructure[:,k+2]=replace(netstructure[:,k+2],missing=>0)
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
#  使用leftjoin对netstructure和partnet进行合并，保留netstructure的所有内容，依据MergeIndex列，生成目标群类的网络一维结构信息
    netstructure=leftjoin(netstructure,partnet,on=:MergeIndex)
#  用replace对原数据进行修改
    netstructure[:,k+2]=replace(netstructure[:,k+2],missing=>0)
  end
#  用sort!在原数据上进行排序，依据EdgeID列进行排序
  sort!(netstructure,[:EdgeID],rev=(false))
#  利用*符合进行string拼接，利用mkdir创建文件夹，总表写入csv文件
  if writemode=="YES"
    wd=pwd()
    netstructurecsv=wd*"/netstructure.csv"
    CSV.write(netstructurecsv,netstructure)
  end
#  对外返还数据
  return(netstructure)
end