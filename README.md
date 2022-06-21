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
* Or you can install this package by gitee:
```
using Pkg
Pkg.add(PackageSpec(url="https://gitee.com/pandalinux/SpeciesToNetworks.jl"))
```
#### Example
* Obtain the basic network properties of different group classes：
```
using RDatasets,SpeciesToNetworks;
dataframe=dataset("datasets","iris");
groupnetinf=Groups2Netinf(dataframe,5;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)
```
* Get the network node label and edge data of the specified group class：
```
using RDatasets,SpeciesToNetworks;
dataframe=dataset("datasets","iris")
labeldf,edgedf,Q,beginqvalue,newcomnumber,begincomnumber=LabelDf(dataframe,5,"setosa";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)
```
* Create a network based on the specified group class：
```
using SpeciesToNetworks,Graphs,DataFrames,RDatasets;
dataframe=dataset("datasets", "iris");
indexdf,edgedf,net=Group2Net(dataframe,5,"setosa";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)
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
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(species1=[1,2],species2=[3,4],species3=[40,50]);
selectdf1=RmPer(dataframe;setper=0.05);
selectdf2=RmPer(dataframe;setper=0.1);
print(dataframe,selectdf1,selectdf2)
```
##### RmZero
* Function description:
```
`RmZero(dataframe::DataFrame)`
Delete all 0 columns.
# Argument
* `dataframe`:A dataframe includes species number information or species abundance information. The column name is the name of each species, and a row means a sample.
# Return
* `selectdf`:Compared with the input dataframe, the returned dataframe does not change the species data, but which columns are all 0 are deleted.
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(species1=[1,2,3],species2=[0,0,0],species3=[-1,0,1]);
selectdf=RmZero(dataframe);
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
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(species1=[1,0,1],species2=[8,5,0],species3=[1,5,9]);
perdf=Per(dataframe);
print(dataframe,perdf)
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
```
* Example:
```
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
```
* Example:
```
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
```
* Example:
```
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
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames;
dataframe=DataFrame(species1=[1,1,0,0,0],species2=[3,3,2,2,2],species3=[1,1,2,2,2],species4=[1,2,3,4,5]);
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg")
```
##### CP2Link
* Function description:
```
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
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05)
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
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Edge2Graph(edgedf,indexdf)
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
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf)
```
##### Avd
* Function description:
```
`Avd(dataframe::DataFrame)`
Calculate average variation degree by the [paper](https://www.researchgate.net/publication/348927379_Specialized_metabolic_functions_of_keystone_taxa_sustain_soil_microbiome_stability).
# Argument
* `dataframe`:A dataframe that stores species abundance data. The column name is the name of each species, and a row means a sample. 
# Return
* `avdvalue`:The value of average variation degree.
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
avdvalue=Avd(dataframe)
```
##### WAvd
* Function description:
```
`WAvd(dataframe::DataFrame)`
Calculate weighted average variation degree.
# Argument
* `dataframe`:A dataframe that stores species abundance data. The column name is the name of each species, and a row means a sample. 
# Return
* `wavd`:The value of weighted average variation degree.
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
wavd=WAvd(dataframe)
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
* `netinf`:A dataframe includes basic network properties.
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf);
netinf=NetInf(net,edgedf,100)
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
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf);
netinfvalue=NetInfValue(net,edgedf,100)
```
##### Group2Net
* Function description:
```
`Group2Net(dataframe::DataFrame,groupcol::Int,groupname::String;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)`
Create a network based on one group.
# Argument
* `dataframe`:A dataframe containing species abundance information and the group of each sample.
* `groupcol`:The index of the column about groups information.
* `groupname`:A group name you want to study.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
* `adjustment`:Select a method to adjust p value, you can use "Bonferroni","BenjaminiHochberg","BenjaminiYekutieli","BenjaminiLiu","Hochberg","Holm","Hommel","Sidak","ForwardStop","BarberCandes","raw".
* `abscorrelation`:Set the judgment conditions of edge connection and require the minimum absolute value of correlation coefficient.
* `pvalue`:Set the judgment conditions of edge connection and require the minimum value of p value.
# Return
* `indexdf`:Generate index numbers for species names.
* `edgedf`:A dataframe for storing edge information.
* `net`:A network based on Graphs.jl.
```
* Example:
```
# Example
using SpeciesToNetworks,Graphs,DataFrames,RDatasets;
dataframe=dataset("datasets", "iris");
indexdf,edgedf,net=Group2Net(dataframe,5,"setosa";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)
```
##### Groups2Netinf
* Function description:
```
`Groups2Netinf(dataframe::DataFrame,groupcol::Int;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,labeln=100)`
Quickly obtain the basic network information of different groups according to the species abundance dataframe.
# Argument
* `dataframe`:A dataframe containing species abundance information and the group of each sample.
* `groupcol`:The index of the column about groups information.
* `method`:You can choose "spearman" or "pearson", these are two algorithms for correlation calculation.
* `adjustment`:Select a method to adjust p value, you can use "Bonferroni","BenjaminiHochberg","BenjaminiYekutieli","BenjaminiLiu","Hochberg","Holm","Hommel","Sidak","ForwardStop","BarberCandes","raw".
* `abscorrelation`:Set the judgment conditions of edge connection and require the minimum absolute value of correlation coefficient.
* `pvalue`:Set the judgment conditions of edge connection and require the minimum value of p value.
# Return
* `groupnetinf`:A dataframe includes basic network properties with different groups.
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("datasets", "iris");
groupnetinf=Groups2Netinf(dataframe,5;method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05,labeln=100)
```
##### ReduceLabel
* Function description:
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
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets,Graphs;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf);
comid=Graphs.label_propagation(net)[1];
newcomid,Q,qvalue,newcomnumber,comnumber=ReduceLabel(comid,net)
```
##### MinLabel
* Function description:
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
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets,Graphs;
dataframe=dataset("MASS","Boston");
linkcor,linkp=SpeciesCP(dataframe,"spearman","BenjaminiHochberg");
indexdf,edgedf,idnetmatrixdf,idnetbooldf=CP2Link(linkcor,linkp,0.6,0.05);
net=Bool2Graph(idnetbooldf);
comid=Graphs.label_propagation(net)[1];
newcomid,Q,beginqvalue,newcomnumber,begincomnumber=MinLabel(comid,net)
```
##### LabelDf
* Function description:
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
```
* Example:
```
# Example
using SpeciesToNetworks,DataFrames,RDatasets;
dataframe=dataset("datasets", "iris");
labeldf,edgedf,Q,beginqvalue,newcomnumber,begincomnumber=LabelDf(dataframe,5,"setosa";method="spearman",adjustment="BenjaminiHochberg",abscorrelation=0.6,pvalue=0.05)
```
