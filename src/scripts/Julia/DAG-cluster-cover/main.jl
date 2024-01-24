include("MyGraph.jl")
include("CVRPInstanceModule.jl")
include("DagFlowInstance.jl")
include("CVRPGraph.jl")
include("BPPGraph.jl")
include("DAGIntFlow.jl")
include("DAGFlowColGen.jl")



using ..CVRPInstanceModule
using ..DagFlowInstanceModule
using ..DAGIntFlow
using ..DAGFlowColGen
using ..CVRPGraph
using ..BPPGraph

function mainCVRPIP()
    instance = CVRPInstanceModule.readInstance(
        "Data/CVRP/Data/A-n32-k5.vrp",
        true,
    )
    #instance = CVRPInstanceModule.readInstance()
    dagFlowInstance = CVRPGraph.constructCVRPGraph(instance)
    DAGIntFlow.setupIP(dagFlowInstance)
end

function mainCVRPColGen(filename = "Data/CVRP/Data/A-n32-k5.vrp")
    instance = CVRPInstanceModule.readInstance(filename, true)
    #instance = CVRPInstanceModule.readInstance("Data/CVRP/Data/P-n19-k2.vrp", true)
    #instance = CVRPInstanceModule.readInstance("Data/CVRP/rand7_1.txt", true)
    #instance = CVRPInstanceModule.readInstance("Data/CVRP/Data/E-n101-k14.vrp", true)
    #instance = CVRPInstanceModule.readInstance("Data/CVRP/Data/E-n101-k8.vrp", true)
    dagFlowInstance = CVRPGraph.constructCVRPGraph(instance)
    DAGFlowColGen.colgen(dagFlowInstance)
end


function mainBPP()
    instance = BPPGraph.readBPP("Data/BPP/Schwerin/Schwerin_1/Schwerin1_BPP1.txt")
    dagFlowInstance = BPPGraph.constructBPPGraph(instance)
    DAGIntFlow.setupIP(dagFlowInstance)
end

function mainBPPColGen()
    #instance = BPPGraph.readBPP("Data/BPP/Irnich_BPP/csAB125_1.txt")
    instance = BPPGraph.readBPP("Data/BPP/Schwerin/Schwerin_1/Schwerin1_BPP1.txt")
    dagFlowInstance = BPPGraph.constructBPPGraph(instance)
    DAGFlowColGen.colgen(dagFlowInstance)
end

function mainGenGraph(FileName::String)
    dagFlowInstance = DagFlowInstanceModule.readDAGFlowArcsInstance(FileName)
    DAGIntFlow.setupIP(dagFlowInstance, FileName)
end

mainGenGraph("Data/wt007_graph_bis/wt007_081_2.dat")
