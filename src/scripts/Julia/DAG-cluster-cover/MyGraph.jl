module MyGraphModule

export MyGraph, MyGraphOut, addNode!, addArcNoDup!, allArcs, outArcs, inArcs, shortestPath, getArc

mutable struct MyArc
    from::Int32
    to::Int32
    cost::Float64
end

mutable struct MyGraph
    # outnodes[i] are the nodes that can be reached using outgoing arcs from node i
    # outCost[i][j] is the cost of the arc (i,outnodes[i][j])
    outArcs::Vector{Vector{Int32}}
    # inNodes[i] are the nodes that can be reached using incoming arcs to node i
    # inCost[i][j] is the cost of the arc (inNodes[i][j],i)
    inArcs::Vector{Vector{Int32}}
    arcs::Vector{MyArc}
    arcDict::Dict{Pair{Int32,Int32},Int64}
end
MyGraph() = MyGraph(Vector{Vector{Int32}}(), Vector{Vector{Int32}}(), Vector{MyArc}(), Dict{Pair{Int32,Int32},Int64}())
MyGraph(n::Int) = MyGraph([Vector{Int32}() for i = 1:n], [Vector{Int32}() for i = 1:n], Vector{MyArc}(), Dict{Pair{Int32,Int32},Int64}())

# returns the pair (|V|,|A|)
function size(G::MyGraph)
    return length(G.outArcs), length(G.arcs)
end

function getArc(g, arcId)
    return g.arcs[arcId]
end

function addNode!(G::MyGraph)
    push!(G.outArcs, Vector{Int32}())
    push!(G.inArcs, Vector{Int32}())
    # return the index of the new node
    return convert(Int32, length(G.outArcs))
end

# adds an arc. Does not check for duplication
function addArcNoDup!(G::MyGraph,from::Int32,to::Int32, cost::Float64)
    push!(G.arcs, MyArc(from,to,cost))
    push!(G.outArcs[from], length(G.arcs))
    push!(G.inArcs[to], length(G.arcs))
    G.arcDict[Pair(from,to)] = length(G.arcs)
end

function cost(G::MyGraph,from::Int32,to::Int32)
    if haskey(G.arcDict, Pair(from,to))
        return G.arcs[G.arcDict[Pair(from,to)]].cost
    else
        throw("cost(G,from,to): Arc not found")
    end
end

function allArcs(G::MyGraph)
    return G.arcs
end

mutable struct ArcIterator
    g::MyGraph
    arcIds::Vector{Int32}
end

function outArcs(g::MyGraph, i::Int32)
    if i > length(g.outArcs)
        throw("Node is out of range")
    end
    return ArcIterator(g, g.outArcs[i])
end

function inArcs(g::MyGraph, i::Int32)
    if i > length(g.inArcs)
        throw("Node is out of range")
    end
    return ArcIterator(g, g.inArcs[i])
end

function Base.iterate(ai::ArcIterator)
    if length(ai.arcIds) == 0
        return nothing
    else
        arc = ai.g.arcs[ai.arcIds[1]]
    end
    return arc,2
end

function Base.iterate(ai::ArcIterator, nextIndex)
    if nextIndex > length(ai.arcIds)
        return nothing
    else
        arc = ai.g.arcs[ai.arcIds[nextIndex]]
        return arc, nextIndex+1
    end
end

mutable struct PromisingNode
    prev::Int32
    cost::Float64
end

function getPath(node, prevNode)
    path = Vector{Int32}()
    push!(path, node)
    while (prevNode[node] != 0)
        node = prevNode[node]
        push!(path, node)
    end
    return reverse(path)
end

function shortestPath(g::MyGraph, nodeOrder::Vector{Int32}, nodeBonus::Vector{Float64}, source::Int32, sink::Int32, maxNPaths::Int64, eps::Float64)
    n,m = size(g)
    cost = typemax(Float64)*ones(Float64,n)
    prevNode = zeros(Int32,n)
    cost[source] = -nodeBonus[source]
    for i in nodeOrder
        if i != source && i!= sink
            nodeCost = typemax(Float64)
            prev = 0
            for arc in inArcs(g,i)
                if cost[arc.from] + arc.cost < nodeCost
                    nodeCost = cost[arc.from] + arc.cost
                    prev = arc.from
                end
            end
            cost[i] = nodeCost-nodeBonus[i]
            prevNode[i] = prev
            #println("$i: best cost: $(cost[i]), prev: $prev")
        end
    end
    # we've reached the sink node. Here we may want to return multiple paths with negative cost
    bestPaths = Vector{PromisingNode}()
    for arc in inArcs(g,sink)
        if cost[arc.from] + arc.cost - nodeBonus[sink] < eps
            push!(bestPaths, PromisingNode(arc.from, cost[arc.from] + arc.cost - nodeBonus[sink]))
        end
    end
    sort!(bestPaths, by = x -> x.cost)
    index = 1
    paths = Vector{Vector{Int32}}()
    costs = Vector{Float64}()
    while index <= maxNPaths && index <= length(bestPaths)
        push!(paths, getPath(bestPaths[index].prev,prevNode))
        push!(paths[index], sink)
        push!(costs, bestPaths[index].cost)
        index += 1
    end
    #println("paths returned: ", paths)
    #println("costs returned: ", costs)
    return paths, costs
end

###############################################################################
###############################################################################
###############################################################################

mutable struct MyGraphOut
    # outnodes[i] is a vector containing the nodes that can be reached using outgoing arcs from node i
    # outCost[i][j] is the cost of the arc (i,outnodes[i][j])
    outNodes::Vector{Vector{Int32}}
    outCost::Vector{Vector{Float64}}
end

#MyGraphOut() = MyGraphOut(Vector{Vector{Int32}}(), Vector{Vector{Float64}())
MyGraphOut() = MyGraphOut([],[])

# helper for iterator functions. Finds the next arc in our ordering
function nextArc(G::MyGraphOut,(fromNode,index)::Tuple{Int64,Int64})
    while (fromNode <= length(G.outNodes)) && (index > length(G.outNodes[fromNode]))
        index = 1
        fromNode += 1
    end
    return (fromNode,index)
end

# initialize the iterator
# Notice, this function has to "belong to Base", otherwise
# it won't work in a "for a in graph" statement
function Base.iterate(G::MyGraphOut)
    (fromNode,index) = nextArc(G, (1,1))
    if fromNode > length(G.outNodes)
        return nothing
    else
        arc = MyArc(fromNode, G.outNodes[fromNode][index], G.outCost[fromNode][index])
        nextState = nextArc(G, (fromNode,index+1))
        return arc, nextState
    end
end

# iteration step
# Notice, this function has to "belong to Base", otherwise
# it won't work in a "for a in graph" statement
function Base.iterate(G::MyGraphOut, (fromNode,index)::Tuple{Int64,Int64})
    if fromNode > length(G.outNodes)
        return nothing
    else
        arc = MyArc(fromNode, G.outNodes[fromNode][index], G.outCost[fromNode][index])
        nextState = nextArc(G, (fromNode,index+1))
        return arc, nextState
    end
end

function addNode!(G::MyGraphOut)
    push!(G.outNodes, Vector{Int32}())
    push!(G.outCost, Vector{Float64}())
    return convert(Int32, length(G.outNodes))
end

# adds an arc. Does not check for duplication
function addArcNoDup!(G::MyGraphOut,from::Int32,to::Int32, cost::Float64)
    push!(G.outNodes[from], to)
    push!(G.outCost[from], cost)
end

# returns the pair (|V|,|A|)
function size(G::MyGraphOut)
    nArcs = 0
    for i=1:length(G.outNodes)
        nArcs += length(G.outNodes[i])
    end
    return length(G.outNodes), nArcs
end

function getNNodes(G::MyGraphOut)
    return length(G.outNodes)
end

function createMyGraph(g::MyGraphOut)
    g2 = MyGraph()
    n = length(g.outNodes)
    g2.outArcs = [Vector{Int32}() for i=1:n]
    g2.inArcs = [Vector{Int32}() for i=1:n]
    for i=1:n
        for j = 1:length(g.outNodes[i])
            otherNode = g.outNodes[i][j]
            cost = g.outCost[i][j]
            addArcNoDup!(g2,Int32(i),otherNode, cost)
        end
    end
    return g2
end


function testMyGraphOut()
    G = MyGraphOut()
    for i=1:10
        addNode!(G)
    end
    addArcNoDup!(G,Int32(2),Int32(3),3.14)
    addArcNoDup!(G,Int32(2),Int32(5),42.42)
    addArcNoDup!(G,Int32(7),Int32(2),12.34)
    addArcNoDup!(G,Int32(9),Int32(1),12.345)
    #println(iterate(G))
    for arc in G
        println(arc)
    end

    G2 = createMyGraph(G)
    println("number of nodes and arcs in G2: ", size(G2))
    println("Arcs leaving node 2:")
    for arc in outArcs(G2,Int32(2))
        println(arc)
    end
    println("Arcs entering node 1:")
    for arc in inArcs(G2,Int32(1))
        println(arc)
    end

end

function testShortestPath()
    G = MyGraph()
    for i=1:10
        addNode!(G)
    end
    addArcNoDup!(G,Int32(1),Int32(3),3.14)
    addArcNoDup!(G,Int32(1),Int32(10),42.42)
    addArcNoDup!(G,Int32(1),Int32(2),12.34)
    addArcNoDup!(G,Int32(2),Int32(5),12.0)
    addArcNoDup!(G,Int32(3),Int32(5),12.0)
    addArcNoDup!(G,Int32(2),Int32(10),3.0)
    addArcNoDup!(G,Int32(5),Int32(10),1.0)

    nodeBonus = zeros(Float64,10)
    nodeBonus[5] = 2
    #shortestPath(g::MyGraph, nodeOrder::Vector{Int32}, nodeBonus::Vector{Float64}, source::Int32, sink::Int32, maxNPaths::Int64, eps::Float64)
    shortestPath(G, collect(Int32(1):Int32(10)), nodeBonus, Int32(1), Int32(10), 5, 17.0) #-0.00001)

end

#MyGraphModule.testMyGraphOut()
#MyGraphModule.testShortestPath()
end
