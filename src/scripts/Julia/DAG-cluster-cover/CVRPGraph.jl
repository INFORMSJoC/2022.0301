#include("CVRPInstanceModule.jl")
#include("MyGraph.jl")
#include("DagFlowInstance.jl")

module CVRPGraph

using ..CVRPInstanceModule
using ..MyGraphModule
using ..DagFlowInstanceModule

function constructCVRPGraph(instance::CVRPInstance)
    #println(instance.c)
    n = instance.n
    nn = instance.nn
    Q = instance.Q
    demand = instance.demand
    println("n=$n, Q=$Q")
    # nodeDemandMap[i,d] should contain the id of the node in the extended graph corresponding to visiting node i in the original
    # graph on a route with total demand d
    nodeDemandMap = zeros(Int32, n, Q)
    for i=1:n
        if demand[i] <= 0
            println("Node $i has zero or negative demand. We cannot construct extended CVRP graph")
            throw("Cannot construct extended CVRP graph. Demand too low")
        end
        if demand[i] > Q
            println("Node $i has a demand that is higher than the vehicle capacity. We cannot construct extended CVRP graph")
            throw("Cannot construct extended CVRP graph. Demand too high")
        end
    end
    graph=MyGraphOut()
    startNode = addNode!(graph)
    endNode = addNode!(graph)
    for i=1:n
        node = recursiveConstruct(convert(Int32, i), demand[i], graph, nodeDemandMap, instance, endNode)
        addArcNoDup!(graph,startNode,node, instance.c[nn,i])
    end
    nNodes, nArcs = MyGraphModule.size(graph)
    println("Extended graph contains $nNodes nodes and $nArcs arcs.")

    nodeInfoDict = Dict{Int32,String}()
    nodePartition = []
    node2Partion = zeros(Float64,nNodes)
    for i=1:n
        nodeSet = Vector{Int32}()
        partitionId = length(nodePartition)+1
        for q=1:Q
            nodeId = nodeDemandMap[i,q]
            if nodeId != 0
                push!(nodeSet, nodeId)
                push!(nodeInfoDict, nodeId => "node $i acc demand $q")
                node2Partion[nodeId] = partitionId
            end
        end
        push!(nodePartition, nodeSet)
    end
    push!(nodeInfoDict, startNode => "Start")
    push!(nodeInfoDict, endNode => "End")

    # Create DagFlowInstance
    gIO = MyGraphModule.createMyGraph(graph)
    # right now the CVRP instance does not contain a max number of vehicles so we just set it to n.
    maxVeh = n

    nodeOrder=Vector{Int32}()
    push!(nodeOrder,startNode)
    for q=1:Q
        for i=1:n
            if nodeDemandMap[i,q] != 0
                push!(nodeOrder, nodeDemandMap[i,q])
            end
        end
    end

    dfi = DagFlowInstance(gIO, nodePartition, node2Partion, startNode, endNode, nodeOrder, [OutFlowLimit(startNode, maxVeh)], Set(startNode), nodeInfoDict)
    return dfi
end

function recursiveConstruct(node::Int32, totalDemand::Int32, graph::MyGraphOut, nodeDemandMap, instance::CVRPInstance, endNode::Int32)
    n = instance.n
    nn = instance.nn
    Q = instance.Q
    demand = instance.demand
    newNode = addNode!(graph)
    nodeDemandMap[node,totalDemand] = newNode
    # connect new node to the end node:
    addArcNoDup!(graph, newNode, endNode, instance.c[node,nn])
    for i=1:n
        if i != node
            newDemand = totalDemand + demand[i]
            if newDemand <= Q
                otherNode = nodeDemandMap[i,newDemand]
                # have we already constructed this extended node?
                if otherNode == 0
                    # no, we need to do that
                    otherNode = recursiveConstruct(convert(Int32, i),newDemand, graph, nodeDemandMap, instance, endNode)
                end
                # connect extended nodes
                addArcNoDup!(graph, newNode, otherNode, instance.c[node,i])
            end
        end
    end
    return newNode
end

#instance = CVRPInstanceModule.readInstance("Data/CVRP/Data/A-n32-k5.vrp", true)
#instance = CVRPInstanceModule.readInstance("Data/CVRP/Data/M-n200-k16.vrp", true)
#graph = constructCVRPGraph(instance)


end
