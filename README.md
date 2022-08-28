### SpeciesToNetworks.jl
#### Overview
* SpeciesToNetworks. jl is a tool to convert species abundance data into undirected network, the basic principle of the tool is to  judge whether there is a connection according to the Spearman or Pearson. 
* You can use ```?``` to read the document of functions in Julia(REPL), or read the README.md published in the [SpeciesToNetworks.jl](https://github.com/JiangXingChi/SpeciesToNetworks.jl) repository.
#### Install
* This package can be installed via Pkg:
```
using Pkg
Pkg.add("SpeciesToNetworks")
```
* Or you can install this package by GitHub:
```
using Pkg;
Pkg.add(PackageSpec(url="https://github.com/JiangXingChi/SpeciesToNetworks.jl"))
```
#### Example
* Example 1:
```
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(Groups=["a","a","a","a","b","b","b","b"],species1=[2,2,1,1,0,5,7,2],species2=[0,0,0,0,3,2,2,2],species3=[1,1,2,2,6,8,2,2],species4=[0,2,2,4,9,3,4,5]);
groupnetinf=Groups2Net(dataframe,1;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100,writemode="NO")
```
* Example 2:
```
using SpeciesToNetworks,DataFrames;
a=DataFrame(group=["ck","ck","ck","test","test","test"],a1=[1,1,2,2,2,1],a2=[3,3,4,2,2,2],a3=[1,1,2,5,2,2],a4=[1,2,1,3,4,5]);
b=DataFrame(group=["ck","ck","ck","test","test","test"],b1=[2,2,1,1,1,2],b2=[3,4,2,2,3,1]);
groupnetinf=Groups2Bipartite(a,b,1;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100,writemode="NO")
```
* Example 3:
```
using SpeciesToNetworks,Graphs,DataFrames;
net=erdos_renyi(10, 30);
NullTimes(net,100,"NullOne",99)
```
* Example 4:
```
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(Groups=["a","a","a","a","b","b","b","b"],species1=[2,2,1,1,0,5,7,2],species2=[0,0,0,0,3,2,2,2],species3=[1,1,2,2,6,8,2,2],species4=[0,2,2,4,9,3,4,5]);
groupnetinf=Groups2Net(dataframe,1;abscorrelation=0,pvalue=1,writemode="YES");
netstructure=NetStructure(dataframe,1;writemode="YES")
```
* Example 5:
```
using SpeciesToNetworks,DataFrames;
a=DataFrame(group=["ck","ck","ck","test","test","test"],a1=[1,1,2,2,2,1],a2=[3,3,4,2,1,2],a3=[1,1,2,5,2,2],a4=[1,2,1,3,4,5]);
b=DataFrame(group=["ck","ck","ck","test","test","test"],b1=[2,2,1,1,1,2],b2=[3,4,2,2,3,1]);
groupnetinf=Groups2Bipartite(a,b,1;adjustment="raw",abscorrelation=0,pvalue=1,writemode="YES");
bipartitestructure=BipartiteStructure(a,b,1;writemode="YES")
```
#### Function
##### RmPer
* Function description:
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
##### RmZeroVector
* Function description:
```
`RmZeroVector(dataframe::DataFrame)`
Delete all 0 columns.
# Argument
* `dataframe`:A dataframe includes species number information or species abundance information. The column name is the name of each species, and a row means a sample.
# Return
* `selectdf`:Compared with the input dataframe, the returned dataframe does not change the species data, but which columns are all 0 are deleted.
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(species1=[1,2,3],species2=[0,0,0],species3=[-1,0,1]);
selectdf=RmZeroVector(dataframe);
print(dataframe,selectdf)
```
##### AllNoZero
* Function description:
```
`AllNoZero(dataframe::DataFrame)`
Delete column containing 0 element.
# Argument
* `dataframe`:A dataframe includes species number information or species abundance information. The column name is the name of each species, and a row means a sample.
# Return
* `selectdf`:Compared with the input dataframe, the returned dataframe does not change the species data, but which columns are containing 0 are deleted.
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(species1=[1,2,3],species2=[0,0,0],species3=[-1,0,1]);
selectdf=AllNoZero(dataframe);
print(dataframe,selectdf)
```
##### Per
* Function description:
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
##### Ck1
* Function description:
```
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
```
##### Ck1s
* Function description:
```
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
```
##### GroupsMean
* Function description:
```
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
```
##### SpeciesCor
* Function description:
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
##### SpeciesPvalue
* Function description:
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
##### PvalueAdjustment
* Function description:
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
##### SpeciesCP
* Function description:
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
##### CP2Link
* Function description:
```
`CP2Link(linkcor::DataFrame,linkp::DataFrame,abscorrelation,pvalue)`
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
dataframe=DataFrame(species1=[1,1,0,0,0],species2=[3,3,2,2,2],species3=[1,1,2,2,2],species4=[1,2,3,4,5]);
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05)
```
##### Edge2Graph
* Function description:
```
`Edge2Graph(edgedf::DataFrame,indexdf::DataFrame)`
Generate a network based on a dataframe about edge information and a dataframe about vertex information.
# Argument
* `edgedf`:A dataframe about edge information.
* `indexdf`:A dataframe about vertex information.
# Return
* `net`:A network.
# Example
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(species1=[1,1,0,0,0],species2=[3,3,2,2,2],species3=[1,1,2,2,2],species4=[1,2,3,4,5]);
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Edge2Graph(edgedf,indexdf);
collect(edges(net))
```
##### Bool2Graph
* Function description:
```
`Bool2Graph(idnetbooldf::DataFrame)`
Generate a network based on a dataframe for storing the network matrix, whether there is a connection expressed by 0 and 1.
# Argument
* `idnetbooldf`:A dataframe for storing the network matrix, whether there is a connection expressed by 0 and 1.
# Return
* `net`:A network.
# Example
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(species1=[1,1,0,0,0],species2=[3,3,2,2,2],species3=[1,1,2,2,2],species4=[1,2,3,4,5]);
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net1=Edge2Graph(edgedf,indexdf);
net2=Bool2Graph(idnetbooldf);
net1==net2
```
##### NetInf
* Function description:
```
`NetInf(net,edgedf::DataFrame,labeln::Int)`
Calculate some network properties by a network.
# Argument
* `net`:A network based on Graphs.jl.
* `edgedf`:A dataframe about edge information.
* `labeln`:Set the times of running label propagation algorithm.
# Return
* `netinf`:A dataframe includes some network properties.
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(species1=[1,1,0,0,0],species2=[3,3,2,2,2],species3=[1,1,2,2,2],species4=[1,2,3,4,5]);
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf);
netinf=NetInf(net,edgedf,1000)
```
##### NetInfValue
* Function description:
```
`NetInfValue(net,edgedf::DataFrame,labeln::Int)`
Calculate some network properties by a network.
# Argument
* `net`:A network based on Graphs.jl.
* `edgedf`:A dataframe about edge information.
* `labeln`:Set the times of running label propagation algorithm.
# Return
* `netinfvalue`:A dataframe includes some network properties, but there's only one column of values.
# Example
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(species1=[1,1,0,0,0],species2=[3,3,2,2,2],species3=[1,1,2,2,2],species4=[1,2,3,4,5]);
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf);
netinfvalue=NetInfValue(net,edgedf,1000)
```
##### Group2Net
* Function description:
```
`Group2Net(dataframe::DataFrame,groupcol::Int,groupname;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector")`
Create a network based on one group.
# Argument
* `dataframe`:A dataframe containing species abundance information and sample groups information.
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
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(Groups=["a","a","a","a","b","b","b","b"],species1=[2,2,1,1,0,5,7,2],species2=[0,0,0,0,3,2,2,2],species3=[1,1,2,2,6,8,2,2],species4=[0,2,2,4,9,3,4,5]);
indexdf1,edgedf1,idnetbooldf1,net1=Group2Net(dataframe,1,"a";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector");
indexdf2,edgedf2,idnetbooldf2,net2=Group2Net(dataframe,1,"a";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="AllNoZero");
indexdf3,edgedf3,idnetbooldf3,net3=Group2Net(dataframe,1,"a";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="raw");
print(idnetbooldf1,idnetbooldf2,idnetbooldf3)
```
##### Groups2Net
* Function description:
```
`Groups2Net(dataframe::DataFrame,groupcol::Int;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100)`
Quickly obtain the basic network information of different groups according to the species abundance dataframe.
# Argument
* `dataframe`:A dataframe containing species abundance information and the group of each sample.
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
using SpeciesToNetworks,DataFrames,Graphs;
dataframe=DataFrame(Groups=["a","a","a","a","b","b","b","b"],species1=[2,2,1,1,0,5,7,2],species2=[0,0,0,0,3,2,2,2],species3=[1,1,2,2,6,8,2,2],species4=[0,2,2,4,9,3,4,5]);
groupnetinf=Groups2Net(dataframe,1;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,colfun="RmZeroVector",labeln=100,writemode="YES")
```
##### BipartiteCP
* Function description:
```
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
```
##### Group2Bipartite
* Function description:
```
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
```
##### Groups2Bipartite
* Function description:
```
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
```
##### NullZero
* Function description:
```
`NullZero(net)`
Generate 0-order null model from an original network.
# Argument
* `net`:A network based on Graphs.jl.
# Return
* `null0net`:A 0-order null model network.
# Example
using Species2Networks,Graphs;
net=erdos_renyi(10, 30);
print(ne(net)==ne(NullZero(net)));
print(degree(net)==degree(NullZero(net)));
print(collect(edges(net))==collect(edges(NullZero(net))))
```
##### NullOne
* Function description:
```
`NullOne(net)`
Generate 1-order null model from an original network.
# Argument
* `net`:A network based on Graphs.jl.
# Return
* `null1net`:A 1-order null model network.
# Example
using Species2Networks,Graphs;
net=erdos_renyi(10, 30);
print(ne(net)==ne(NullOne(net)));
print(degree(net)==degree(NullOne(net)));
print(collect(edges(net))==collect(edges(NullOne(net))))
```
##### NullInf
* Function description:
```
`NullInf(rawnet,labeln::Int,nullmode::String)`
Generate a corresponding null model network according to an original network, and obtain some properties of the null model network.
# Argument
* `rawnet`:An original network.
* `labeln`:Set the times of running label propagation algorithm.
* `nullmode`:You can choose "NullZero" or "NullOne".
# Return
* `netinf`:A dataframe includes some network properties with a null model network.
# Example
using Species2Networks,Graphs,DataFrames;
net=erdos_renyi(10, 30);
NullInf(net,100,"NullOne")
```
##### NullInfValue
* Function description:
```
`NullInfValue(rawnet,labeln::Int,nullmode::String)`
Generate a corresponding null model network according to an original network, and obtain some properties of the null model network.
# Argument
* `rawnet`:An original network.
* `labeln`:Set the times of running label propagation algorithm.
* `nullmode`:You can choose "NullZero" or "NullOne".
# Return
* `netinfvalue`:A dataframe includes some network properties with a null model network,but there is only value information.
# Example
using Species2Networks,Graphs,DataFrames;
net=erdos_renyi(10, 30);
NullInfValue(net,100,"NullOne")
```
##### NullTimes
* Function description:
```
`NullTimes(rawnet,labeln::Int,nullmode::String,timen::Int)`
Generate some corresponding null model networks according to an original network, and obtain some properties of the null model networks.
# Argument
* `rawnet`:An original network.
* `labeln`:Set the times of running label propagation algorithm.
* `nullmode`:You can choose "NullZero" or "NullOne".
* `timen`:Set the times of running null model algorithm.
# Return
* `nulltimesdf`:A dataframe includes some network properties with some null model networks.
# Example
using SpeciesToNetworks,Graphs,DataFrames;
net=erdos_renyi(10, 30);
NullTimes(net,100,"NullOne",99)
```
##### NetStructure
* Function description:
```
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
```
##### BipartiteStructure
* Function description:
```
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
```