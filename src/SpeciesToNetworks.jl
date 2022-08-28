module SpeciesToNetworks

#  加载需要的包
using DataFrames,MultipleTesting,HypothesisTests,StatsBase,Statistics,Graphs,CategoricalArrays,CSV

#  加载脚本
include("1-datapre.jl")
include("2-datagenerate.jl")
include("3-corcal.jl")
include("4-cornet.jl")
include("5-netinf.jl")
include("6-groupnet.jl")
include("7-bipartite.jl")
include("8-nullbase.jl")
include("9-netstructure.jl")
include("10-bipartitestructure.jl")

#  输出函数
export RmPer,RmZeroVector,AllNoZero,Per,
       Ck1,Ck1s,GroupsMean,
       SpeciesCor,SpeciesPvalue,PvalueAdjustment,SpeciesCP,
       CP2Link,Edge2Graph,Bool2Graph,
       NetInf,NetInfValue,
       Group2Net,Groups2Net,
       BipartiteCP,Group2Bipartite,Groups2Bipartite,
       NullZero,NullOne,NullInf,NullInfValue,NullTimes,
       NetStructure,
       BipartiteStructure

#  设置文档
"""
SpeciesToNetworks. jl is a tool to convert species abundance data into undirected network, the basic principle of the tool is to  judge whether there is a connection according to the Spearman or Pearson. 
"""
SpeciesToNetworks

end
