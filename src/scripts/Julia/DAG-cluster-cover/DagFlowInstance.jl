# include("MyGraph.jl")

module DagFlowInstanceModule

using ..MyGraphModule

export OutFlowLimit, DagFlowInstance, DagFlowInstanceArcs

mutable struct OutFlowLimit
    node::Int32
    limit::Int64
end

mutable struct DagFlowInstance
    g::MyGraph
    nodePartition::Vector{Vector{Int32}}
    #node2Partition[i] gives the partition of node i
    node2Partion::Vector{Int32}
    sourceNode::Int32
    sinkNode::Int32
    nodeOrder::Vector{Int32}
    outFlowLimits::Vector{OutFlowLimit}
    setOutFlowNodes::Set{Int32}
    nodeInfoDict::Dict{Int32,String}
end

mutable struct DagFlowInstanceArcs
    g::MyGraph
    arcPartition::Vector{Vector{Int64}}
    #arc2Partition[i] gives the partition of arc i
    arc2Partion::Vector{Int32}
    sourceNode::Int32
    sinkNode::Int32
    nodeOrder::Vector{Int32}
    outFlowLimits::Vector{OutFlowLimit}
    setOutFlowNodes::Set{Int32}
    nodeInfoDict::Dict{Int32,String}
    highArcs::Vector{Int64}
    lowArcs::Vector{Int64}
end

function readDAGFlowArcsInstance(fileName::String)
    open(fileName) do f
      allFile = read(f, String)
      allNumbers = map(x -> parse(Int64,x), split(allFile))
      #println(allNumbers)
      nNodes = allNumbers[1]
      nArcs = allNumbers[2]
      nPartitions = allNumbers[3]
      nMaxPaths = allNumbers[4]
      graph=MyGraph(nNodes)
      index=5
      for i=1:nArcs
          src = allNumbers[index]   +1
          dst = allNumbers[index+1] +1
          cost = allNumbers[index+2]
          #println("src: $src, dst: $dst, cost: $cost")
          addArcNoDup!(graph, Int32(src), Int32(dst), Float64(cost))
          index += 3
      end
      arc2Partition = zeros(Int64,nArcs)
      arcPartition = Vector{Vector{Int64}}()
      highArcsSet = BitSet()
      for i=1:nPartitions
          nArcsInPartition = allNumbers[index]
          #println("nArcsInPartition: $nArcsInPartition")
          index += 1
          arcs = Vector{Int64}()
          for j=1:nArcsInPartition
              arcId = allNumbers[index]+1
              push!(arcs,arcId)
              arc2Partition[arcId] = i
              push!(highArcsSet, arcId)
              index+=1
          end
          #println("arcs $i: ", arcs)
          push!(arcPartition,arcs)
      end
      nodeInfoDict = Dict{Int32,String}()
      for i=1:nNodes
          node = allNumbers[index]
          time = allNumbers[index+1]
          #println("$i: $node $time")
          index += 2
          push!(nodeInfoDict, i => "$node $time")
      end

      lowArcs = Vector{Int64}()
      for arcId=1:nArcs
          if !(arcId in highArcsSet)
              push!(lowArcs, arcId)
          end
      end

      ofl = OutFlowLimit(1,nMaxPaths)
      setOutFlowNodes = Set{Int32}(Int32(1))
      return DagFlowInstanceArcs(graph, arcPartition, arc2Partition, 1, nNodes, [], [ofl], setOutFlowNodes, nodeInfoDict, collect(highArcsSet), lowArcs )
    end
end

#println(readDAGFlowArcsInstance("data/graphs/wt008_028_2.txt"))

end
