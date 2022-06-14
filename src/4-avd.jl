#  加载包
using DataFrames,StatsBase,Statistics

#  定义一个计算AVD值的函数
"""
`Avd(dataframe::DataFrame)`
Calculate average variation degree by the [paper](https://www.researchgate.net/publication/348927379_Specialized_metabolic_functions_of_keystone_taxa_sustain_soil_microbiome_stability).
# Argument
* `dataframe`:A dataframe that stores species abundance data. The column name is the name of each species, and a row means a sample. 
# Return
* `avdvalue`:The value of average variation degree.
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
avdvalue=Avd(dataframe)
"""
function Avd(dataframe::DataFrame)
#  去除不存在的列
  dataframe=RmZero(dataframe)
#  获取dataframe的尺寸
  a,b=size(dataframe)
#  制作一个储存每列的ai和的数据结构
  aicolsum=fill(0.0,b)
#  根据ai_layout的长度进行内部循环
  for i in 1:b
#  计算第i列的标准化数据
    standinf=StatsBase.fit(ZScoreTransform,float.(dataframe[:,i]);center=true,scale=true)
#  计算获得转换后的数据，用点操作对float进行广播
    standcoli=StatsBase.transform(standinf,float.(dataframe[:,i]))
#  利用点操作进行广播语法，套上绝对值
    standabsi=abs.(standcoli)
#  计算转换数据的和
    colaisum=sum(standabsi)
#  把和的结果填充进aicolsum
    aicolsum[i]=colaisum
  end
#  计算AVD值
  avdvalue=sum(aicolsum)/(a*b)
#  对外部返回值
  return(avdvalue)
end

#  定义一个计算AVD值的函数
"""
`WAvd(dataframe::DataFrame)`
Calculate weighted average variation degree.
# Argument
* `dataframe`:A dataframe that stores species abundance data. The column name is the name of each species, and a row means a sample. 
# Return
* `wavd`:The value of weighted average variation degree.
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
wavd=WAvd(dataframe)
"""
function WAvd(dataframe::DataFrame)
#  去除不存在的列
  dataframe=RmZero(dataframe)
#  获取dataframe的尺寸
  a,b=size(dataframe)
#  制作一个储存每列的ai和的数据结构
  aicolsum=fill(0.0,b)
#  根据ai_layout的长度进行内部循环
  for i in 1:b
#  计算第i列的标准化数据
    standinf=StatsBase.fit(ZScoreTransform,float.(dataframe[:,i]);center=true,scale=true)
#  计算获得转换后的数据，用点操作对float进行广播
    standcoli=StatsBase.transform(standinf,float.(dataframe[:,i]))
#  利用点操作进行广播语法，套上绝对值
    standabsi=abs.(standcoli)
#  计算转换数据的和
    colaisum=sum(standabsi)
#  把和的结果填充进aicolsum
    aicolsum[i]=colaisum
  end
#  计算列和，利用点语法进行广播
  colsum=sum.(eachcol(dataframe))
#  把列和再求和，获得总体数目
  allsum=sum(colsum)
#  利用点语法，将列和中的元素一个个除以总体和
  colper=colsum./allsum
#  利用期望的方式，利用点操作进行广播相乘，除以a消除样本数量的影响，此时因为时期望的思想因此不需要除以b
  wavdi=(aicolsum.*colper)/a
#  求wavdi的和
  wavd=sum(wavdi)
#  对外界返回带权重的avd值
  return(wavd)
end