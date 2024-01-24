#include("MyGraph.jl")
#include("DagFlowInstance.jl")

module BPPGraph

using ..MyGraphModule
using ..DagFlowInstanceModule

mutable struct BPPinstance
    capacity::Int64
    demands::Vector{Int64}
end
BPPinstance() = BPPinstance(0,[])

function readBPP(fileName::String)
    instance = BPPinstance()
    open(fileName) do f
        allFile = read(f, String)
        allNumbers = map(x -> parse(Int64,x), split(allFile))
        n = allNumbers[1]
        instance.capacity = allNumbers[2]
        for i=1:n
            push!(instance.demands, allNumbers[2+i])
        end
    end
    return instance
end

mutable struct BPPNode
    nodeId::Int32
    accDemand::Int64
    item::Int32
end

function constructBPPGraph(instance::BPPinstance)
    graph=MyGraph()
    startNode = addNode!(graph)
    endNode = addNode!(graph)

    allBPPNodes = Vector{BPPNode}()
    push!(allBPPNodes, BPPNode(startNode,0,0))
    n = length(instance.demands)
    cap = instance.capacity
    println("#Items: $n")
    println("Capacity: $cap")

    nodePartition = []
    for i=1:n
        println("Adding nodes and arcs for item $i")
        nodeSet = Vector{Int32}()
        demand = instance.demands[i]
        dictDemandToNode = Dict{Int64,Int32}()
        bppNodesForItem = Vector{BPPNode}()
        for bppNode in allBPPNodes
            newAccDemand = bppNode.accDemand + demand
            #println("newAccDemand = $newAccDemand")
            if newAccDemand <= cap
                # Do we already have a node in the graph corresponding to item i and accumulated demand <newAccDemand>
                if !haskey(dictDemandToNode, newAccDemand)
                    nodeId = addNode!(graph)
                    push!(bppNodesForItem, BPPNode(nodeId,newAccDemand, i))
                    push!(dictDemandToNode, newAccDemand => nodeId)
                    # connect this node to the end node
                    addArcNoDup!(graph, nodeId, endNode, 1.0)
                    push!(nodeSet, nodeId)
                else
                    nodeId = get(dictDemandToNode, newAccDemand, -1)
                end
                addArcNoDup!(graph, bppNode.nodeId, nodeId, 0.0)
            end
        end
        # add newly created nodes to <allBPPNodes>
        # (we cannot add them on fly since we would then go back and iterate over them at "for bppNode in allBPPNodes" and then create unwanted "duplicates" )
        for bppNode in bppNodesForItem
            push!(allBPPNodes, bppNode)
        end
        push!(nodePartition, nodeSet)
    end

    nNodes, nArcs = MyGraphModule.size(graph)
    println("Extended graph contains $nNodes nodes and $nArcs arcs.")

    sort(allBPPNodes, by = x -> x.accDemand)
    nodeOrder=Vector{Int32}()
    nodeInfoDict = Dict{Int32,String}()
    node2Partion = zeros(Float64,nNodes)

    push!(nodeOrder,startNode)
    for bppNode in allBPPNodes
        push!(nodeOrder,bppNode.nodeId)
        node2Partion[bppNode.nodeId] = bppNode.item
        push!(nodeInfoDict, bppNode.nodeId => "item $(bppNode.item) acc demand $(bppNode.accDemand)")
    end
    push!(nodeInfoDict, startNode => "Start")
    push!(nodeInfoDict, endNode => "End")

    dfi = DagFlowInstance(graph, nodePartition, node2Partion, startNode, endNode, nodeOrder, [], Set{Int32}(), nodeInfoDict )
    return dfi
end

end
