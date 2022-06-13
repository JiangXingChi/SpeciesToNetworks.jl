#  加载包
using MultipleTesting,HypothesisTests,StatsBase,Statistics,DataFrames

#  设置一个计算关系系数的函数
```
`SpeciesCor(x::Vector,y::Vector,method::String)`
To calculate the correlation coefficient between x vector and y vector, you need to specify whether Pearson method or Spearman method is used.
# Argument
* `x`:x is a vector, its length must be equal to y.
* `y`:y is a vector, its length must be equal to x.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
# Return
* `i`:Return the correlation coefficient of vector x and vector y.
# Example
using SpeciesToNetworks,DataFrames;
x=[1,2,3];
y=[1,2,300];
i1=SpeciesCor(x,y,"spearman");
i2=SpeciesCor(x,y,"pearson");
print(i1,i2)
```
function SpeciesCor(x::Vector,y::Vector,method::String)
#  当目的是进行spearman时，调用目的函数corspearman
  if method=="spearman"
    i=StatsBase.corspearman(x,y)
#  当目的时进行pearson时，调用目的函数cor
  elseif method=="pearson"
    i=Statistics.cor(x,y)
  end
#  对外界返回相关系数
  return(i)
end

#  设置一个计算关系系数的原始p值的函数
```
`SpeciesPvalue(x::Vector,y::Vector,method::String)`
To calculate the correlation p value between x vector and y vector, you need to specify whether Pearson method or Spearman method is used.
# Argument
* `x`:x is a vector, its length must be equal to y.
* `y`:y is a vector, its length must be equal to x.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
# Return
* `i`:Return the correlation p value of vector x and vector y.
# Example
using SpeciesToNetworks,DataFrames;
x=[1,2,3];
y=[1,2,300];
i1=SpeciesPvalue(x,y,"spearman");
i2=SpeciesPvalue(x,y,"pearson");
print(i1,i2)
```
function SpeciesPvalue(x::Vector,y::Vector,method::String)
#  当目的时进行spearman时，调用目的函数pvalue和tiedrank
  if method=="spearman"
    i=HypothesisTests.pvalue(HypothesisTests.CorrelationTest(StatsBase.tiedrank(x),StatsBase.tiedrank(y)))
#  当目的时进行pearson时，调用目的函数pvalue
  elseif method=="pearson"
    i=HypothesisTests.pvalue(HypothesisTests.CorrelationTest(x,y))
  end
#  对外界返回原始p值
  return(i)
end

#  设置一个函数对p值集合进行矫正
#  根据adjustment决定矫正方法，调用MultipleTesting的相关函数以及不作任何处理的raw方法
```
`PvalueAdjustment(x::Vector,adjustment::String)`
Adjust multiple p values.
# Argument
* `x`:x is a vector consisting of multiple p values.
* `adjustment`:Select a method to adjust p value, you can use "Bonferroni","BenjaminiHochberg","BenjaminiYekutieli","BenjaminiLiu","Hochberg","Holm","Hommel","Sidak","ForwardStop","BarberCandes","raw".
# Return
* `is`:Return the adjusted p values.
# Example
using SpeciesToNetworks;
x=[0.05,0.06,0.12,0.07,0.23,0.89,0.43,0.08,0.16];
is1=PvalueAdjustment(x,"raw");
is2=PvalueAdjustment(x,"BenjaminiHochberg");
print(x,is1,is2)
```
function PvalueAdjustment(x::Vector,adjustment::String)
  if adjustment=="Bonferroni"
    is=MultipleTesting.adjust(MultipleTesting.PValues(x),Bonferroni())
  elseif adjustment=="BenjaminiHochberg"
    is=MultipleTesting.adjust(MultipleTesting.PValues(x),BenjaminiHochberg())
  elseif adjustment=="BenjaminiYekutieli"
    is=MultipleTesting.adjust(MultipleTesting.PValues(x),BenjaminiYekutieli())
  elseif adjustment=="BenjaminiLiu"
    is=MultipleTesting.adjust(MultipleTesting.PValues(x),BenjaminiLiu())
  elseif adjustment=="Hochberg"
    is=MultipleTesting.adjust(MultipleTesting.PValues(x),Hochberg())
  elseif adjustment=="Holm"
    is=MultipleTesting.adjust(MultipleTesting.PValues(x),Holm())
  elseif adjustment=="Hommel"
    is=MultipleTesting.adjust(MultipleTesting.PValues(x),Hommel())
  elseif adjustment=="Sidak"
    is=MultipleTesting.adjust(MultipleTesting.PValues(x),Sidak())
  elseif adjustment=="ForwardStop"
    is=MultipleTesting.adjust(MultipleTesting.PValues(x),ForwardStop())
  elseif adjustment=="BarberCandes"
    is=MultipleTesting.adjust(MultipleTesting.PValues(x),BarberCandes())
  elseif adjustment=="raw"
    is=x
  end
#  对外界返回矫正p值的集合
  return(is)
end

#  datafram时iris(只有数字)那样的数据集
#  method支持spearman和pearson
#  adjustment支持Bonferroni,BenjaminiHochberg,
#  BenjaminiHochbergAdaptive,BenjaminiYekutieli,
#  BenjaminiLiu,Hochberg
#  Holm,Hommel,Sidak,ForwardStop
#  BarberCandes
```
`SpeciesCP(dataframe::DataFrame,method::String,adjustment::String)`
The abundance data of species can be transformed into two dataframes to store the correlation coefficient and correlation p value respectively.
# Argument
* `dataframe`:A dataframe that stores species abundance data. The column name is the name of each species, and a row means a sample.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
* `adjustment`:Select a method to adjust p value, you can use "Bonferroni","BenjaminiHochberg","BenjaminiYekutieli","BenjaminiLiu","Hochberg","Holm","Hommel","Sidak","ForwardStop","BarberCandes","raw".
# Return
* `linkcor`:A dataframe that stores the correlation coefficient.
* `linkp`:A dataframe that stores the correlation p value.
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(species1=[1,1,0,0,0],species2=[3,3,2,2,2],species3=[1,1,2,2,2],species4=[1,2,3,4,5]);
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg")
```
function SpeciesCP(dataframe::DataFrame,method::String,adjustment::String)
#  获取dataframe信息，实际上a就是样本数，b就是物种数/特征数
  a,b=size(dataframe)
#  实际需要计算的p值个数，用Int将其转化为整数
  cpnumber=Int(b*(b-1)/2)
#  创建空白数组array，用来储存结果
  corvalue=fill(0.0,(cpnumber))
  pvalue=fill(0.0,(cpnumber))
#  通过循环进行计算cor值和原始p值
  k=1
  for i in 1:(b-1)
    for j in (i+1):b
      corvalue[k]=SpeciesCor(dataframe[:,i],dataframe[:,j],method)
      pvalue[k]=SpeciesPvalue(dataframe[:,i],dataframe[:,j],method)
      k=k+1
    end
  end
#  对原始p值进行矫正
  adjpvalue=PvalueAdjustment(pvalue,adjustment)
#  创建cor数据框与adjp数据框，用以储存结果
#  这里实际上并没计算自相关，所以linkcor储存相关性时初始化值为1.0，实际上就是自相关的相关性值
#  同理，linkp储存矫正p值时初始化值为0.0，因为自相关的p值必然为显著
  linkcor=DataFrame(fill(1.0,(b,b)),propertynames(dataframe))
  linkp=DataFrame(fill(0.0,(b,b)),propertynames(dataframe))
#  通过循环将计算结果写入数据框
  k=1
  for i in 1:(b-1)
    for j in (i+1):b
      linkcor[i,j]=corvalue[k]
      linkcor[j,i]=corvalue[k]
      linkp[i,j]=adjpvalue[k]
      linkp[j,i]=adjpvalue[k]
      k=k+1
    end
  end
#  创建id数据框，并进行数据框合并
  linkcor=hcat(DataFrame(CorID=names(dataframe)),linkcor)
  linkp=hcat(DataFrame(PvalueID=names(dataframe)),linkp)
#  返回计算结果
  return(linkcor,linkp)
end