module SpeciesToNetworks

#  加载需要的包
using DataFrames,MultipleTesting,HypothesisTests,StatsBase,Statistics,Graphs,CategoricalArrays

#  加载脚本
include("1-datapre.jl")
include("2-corcal.jl")
include("3-cornet.jl")
include("4-avd.jl")
include("5-netinf.jl")
include("6-groupnet.jl")
include("7-community.jl")

#  输出函数
export RmPer,RmZero,Per,
       SpeciesCor,SpeciesPvalue,PvalueAdjustment,SpeciesCP,
       CP2Link,Edge2Graph,Bool2Graph,
       Avd,WAvd,
       NetInf,NetInfValue,
       Group2Net,Groups2Netinf,
       ReduceLabel,MinLabel,LabelDf

#  设置文档
"""
SpeciesToNetworks. jl is a tool to convert species abundance data into undirected network, the basic principle of the tool is to  judge whether there is a connection according to the Spearman or Pearson. 
"""
SpeciesToNetworks

end
