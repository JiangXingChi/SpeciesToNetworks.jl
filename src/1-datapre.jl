#  加载包
using DataFrames,Statistics

#  用于删除总丰度不足0.01%的物种
```
`RmPer(dataframe::DataFrame;setper=0.01)`
Delete the species whose total abundance does not reach the set abundance, and the default minimum total abundance is 0.01.
# Argument
* `dataframe`:A dataframe includes species number information or species abundance information. The column name is the name of each species, and a row means a sample.
* `setper`:Set a number to control the minimum percentage,the default value is 0.01,this means species with total abundance less than 1% were deleted.
# Return
* `selectdf`:Compared with the input dataframe, the returned dataframe does not change the species data, but deletes the columns of species whose abundance does not reach the preset minimum total abundance.
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(species1=[1,2],species2=[3,4],species3=[40,50]);
selectdf1=RmPer(dataframe;setper=0.05);
selectdf2=RmPer(dataframe;setper=0.1);
print(dataframe,selectdf1,selectdf2)
```
function RmPer(dataframe::DataFrame;setper=0.01)
#  计算列和，利用点语法进行广播
    colsum=sum.(eachcol(dataframe))
#  把列和再求和，获得总体数目
    allsum=sum(colsum)
#  利用点语法，将列和中的元素一个个除以总体和
    colper=colsum./allsum
#  把总丰度数据利用dataframe进行包装，附上列信息 
    allper=DataFrame(ID=propertynames(dataframe),Percentage=colper)
#  筛选总丰度大于设定值的行，利用点语法一个个的判定
    selectper=allper[allper.Percentage .> setper,:]
#  通过propertyname，把原数据框中符合要求的列筛选出来
    selectdf=dataframe[:,selectper.ID]
#  对外界返回整理数据
    return(selectdf)
  end

#  用以处理各个分群，把独自类群中不存在的物种剔除
```
`RmZero(dataframe::DataFrame)`
Delete all 0 columns.
# Argument
* `dataframe`:A dataframe includes species number information or species abundance information. The column name is the name of each species, and a row means a sample.
# Return
* `selectdf`:Compared with the input dataframe, the returned dataframe does not change the species data, but which columns are all 0 are deleted.
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(species1=[1,2,3],species2=[0,0,0],species3=[-1,0,1]);
selectdf=RmZero(dataframe);
print(dataframe,selectdf)
```
function RmZero(dataframe::DataFrame)
#  获取dataframe的一维二维长度
  a,b=size(dataframe)
#  创建一个与行数长度相同的0.0向量
  zerov=fill(0.0,a)
#  创建一个存储需要删除列索引和列名称的数据框
  delnamedf=DataFrame(Index=Int[],Names=String[])
#  获取数据框的列名
  namev=names(dataframe)
#  从第一列到最后一列极限循环
  for i in 1:b
#  当循环所在列与创建的零向量相等时，把此时的列数和名称推送到delnamedf中
    if dataframe[:,i] == zerov
      push!(delnamedf,(i,namev[i]))
    end
  end
#  用Not删除全部为0向量的列
  newdf=dataframe[:,Not(delnamedf.Index)]
#  对外界返回删除全为0列的数据框
  return(newdf)
end

#  用于获取丰度矩阵
```
`Per(dataframe::DataFrame)`
Convert the dataframe recording the number of species into the species abundance dataframe.
# Argument
* `dataframe`:A dataframe containing only species number information.The column name is the name of each species, and a row means a sample.
# Return
* `perdf`:A dataframe containing only species abundance information. The column name is the name of each species, and a row means a sample.
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(species1=[1,0,1],species2=[8,5,0],species3=[1,5,9]);
perdf=Per(dataframe);
print(dataframe,perdf)
```
function Per(dataframe::DataFrame)
#  计算行和，利用点语法进行广播
  rowsum=sum.(eachrow(dataframe))
#  获取数据尺寸
  a,b=size(dataframe)
#  制造一个数据框，用于储存结果
  perdf=DataFrame(fill(0.0,(a,b)),propertynames(dataframe))
#  通过循环，一列一列的计算
  for i in 1:b
#  利用点语法，将dataframe的目标列的元素一个个除以rowsum中的元素
    perdf[:,i]=dataframe[:,i]./rowsum
  end
#  对外界返回数据
  return(perdf)
end

