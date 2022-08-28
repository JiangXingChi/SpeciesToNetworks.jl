#  加载包
using CSV,DataFrames,CategoricalArrays,Statistics

#  处理单独群类数据，生成新数据
"""
`Ck1(dataframe::DataFrame,groupcol::Int)`
Generate new data by one group.
# Argument
* `dataframe`:A dataframe containing species abundance information and sample groups information.
* `groupcol`:The index of the column about groups information.
# Return
*`newdataframe`:A new dataframe containing species abundance information and sample groups information.
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(groups=["a","a","a","a"],species1=[1,2,3,4],species2=[0,0,0,0],species3=[1,1,2,2]);
newdataframe=Ck1(dataframe,1)
"""
function Ck1(dataframe::DataFrame,groupcol::Int)
#  获取群类信息
  groupinf=levels(CategoricalArray(dataframe[:,groupcol]))
#  获取群类个数
  groupnumber=Int(size(dataframe,1))
#  获取新群类个数
  newgroupnumber=groupnumber-1
#  初始化处理
  i=1
#  复制数据
  tempdataframe=copy(dataframe)
#  原先的数据库删除第i行
  tempdataframe=tempdataframe[Not(i),:]
#  生成新的群类名称
  newgroupinf=groupinf[1]*" delete "*"index"*string(i)
#  根据新的群类名称生成名称向量
  tempgroupname=repeat([newgroupinf],newgroupnumber)
#  把名称向量替换到对应的列
  tempdataframe[:,groupcol]=tempgroupname
#  初始化新数据框
  newdataframe=copy(tempdataframe)
#  进入循环处理
  for i in 2:groupnumber
#  复制数据
    tempdataframe=copy(dataframe)
#  原先的数据库删除第i行
    tempdataframe=tempdataframe[Not(i),:]
#  生成新的群类名称
    newgroupinf=groupinf[1]*" delete "*"index"*string(i)
#  根据新的群类名称生成名称向量
    tempgroupname=repeat([newgroupinf],newgroupnumber)
#  把名称向量替换到对应的列
    tempdataframe[:,groupcol]=tempgroupname
#  合并数据库
    newdataframe=vcat(newdataframe,tempdataframe)
  end
#  对外返回数据
  return(newdataframe)
end

#   处理多群类数据，生成新数据
"""
`Ck1s(dataframe::DataFrame,groupcol::Int)`
Generate new data by groups.
# Argument
* `dataframe`:A dataframe containing species abundance information and sample groups information.
* `groupcol`:The index of the column about groups information.
# Return
*`generatedf`:A new dataframe containing species abundance information and sample groups information.
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(Groups=["a","a","a","a","b","b","b","b"],species1=[2,2,1,1,0,5,7,2],species2=[0,0,0,0,3,2,2,2],species3=[1,1,2,2,6,8,2,2],species4=[0,2,2,4,9,3,4,5]);
generatedf=Ck1s(dataframe,1)
"""
function Ck1s(dataframe::DataFrame,groupcol::Int)
#  获取群类信息
  groupinf=levels(CategoricalArray(dataframe[:,groupcol]))
#  确认有几个群
  imax=length(groupinf)
#  初始化循环
  i=1
#  筛选满足某个种群名称的行
  tempdf=dataframe[dataframe[:,groupcol].==groupinf[i],:]
#  利用Ck1函数生成数据
  generatedf=Ck1(tempdf,groupcol)
#  开始启动循环
  for i in 2:imax
#  筛选满足某个种群名称的行
    tempdf=dataframe[dataframe[:,groupcol].==groupinf[i],:]
#  利用Ck1函数生成数据
    tempgeneratedf=Ck1(tempdf,groupcol)
#  合并数据
    generatedf=vcat(generatedf,tempgeneratedf)
  end
#  对外返回数据
  return(generatedf)
end

#  获取物种平均数据
"""
`GroupsMean(dataframe::DataFrame,groupcol::Int)`
Generate species average data.
# Argument
* `dataframe`:A dataframe containing species abundance information and sample groups information.
* `groupcol`:The index of the column about groups information.
# Return
*`meandf`:A new dataframe containing average species abundance information.
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(Groups=["a","a","a","a","b","b","b","b"],species1=[2,2,1,1,0,5,7,2],species2=[0,0,0,0,3,2,2,2],species3=[1,1,2,2,6,8,2,2],species4=[0,2,2,4,9,3,4,5]);
meandf=GroupsMean(dataframe,1)
"""
function GroupsMean(dataframe::DataFrame,groupcol::Int)
#  获取群类信息
  groupinf=levels(CategoricalArray(dataframe[:,groupcol]))
#  确认有几个物种
  imax=size(dataframe,2)-1
#  获取列名
  propertynamesvector=propertynames(dataframe)[Not(groupcol)]
#  获取GroupedDataFrame
  groupdf=groupby(dataframe,groupcol)
#  进行循环初始化，从第1个物种开始
  i=1
#  利用combine获取第1个物种的平均值
  meandf=combine(groupdf,(propertynamesvector[i])=>mean)
#  开始进行循环
  for i in 2:imax
#  利用combine获取第1个物种的平均值
    tempmeandf=combine(groupdf,(propertynamesvector[i])=>mean)
#  进行数据合并
    meandf=leftjoin(meandf,tempmeandf,on=propertynames(dataframe)[groupcol])
  end
#  对外返回数据
  return(meandf)
end