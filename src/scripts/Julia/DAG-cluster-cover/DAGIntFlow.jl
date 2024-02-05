# include("MyGraph.jl")
include("DagFlowInstance.jl")

include("../ZeroHalfCuts/ZeroHalfCuts.jl")

module DAGIntFlow

using JuMP
using CPLEX
using Printf
using MathProgBase
using ..MyGraphModule
using ..DagFlowInstanceModule
using ..ZeroHalfChvatGom


#=
function createVarDict(graph::MyGraphOut)
    varDict = Dict{Tuple{Int32,Int32},Int64}()
    n = getNNodes(graph)
    id = 1
    for i=1:n
        for node in graph.outNodes[i]
            varDict[(i,node)] = id
            id += 1
        end
    end
    return varDict
end
=#

function setupIP(instance::DagFlowInstance)
    #varDict = createVarDict(g)
    g = instance.g
    n,m = MyGraphModule.size(g)
    # the nodes where flow conservation is active includes all ndoes except the source and the sink:
    flowConsNodes = setdiff(1:n, [instance.sourceNode instance.sinkNode])

    IP = Model(solver=CplexSolver())

    #gInOut = createMyGraph(graph)
    # @@@ low arcs should be general integers!
    @variable(IP, x[allArcs(g)], Bin)

    @objective(IP, Min, sum(a.cost*x[a] for a in allArcs(g)))

    @constraint(IP, [subset in instance.nodePartition], sum(x[a] for node in subset for a in outArcs(g,node)) == 1)
    #=for subset in instance.nodePartition
        println("subset = $subset")
        for node in subset
            println("Node = $node")
            for a in outArcs(g,node)
                println(a)
            end
            throw("stopping!!!")
        end
    end=#
    @constraint(IP, [node in flowConsNodes], sum(x[a] for a in outArcs(g,Int32(node))) == sum(x[a] for a in inArcs(g,Int32(node))))
    @constraint(IP, [outFlowLimit in instance.outFlowLimits], sum(x[a] for a in outArcs(g,outFlowLimit.node)) <= outFlowLimit.limit)
    writeLP(IP, "DagIntFlow.lp", genericnames=true);
    solve(IP)
    println("objective value: ", getobjectivevalue(IP))
end

function writeDOTFile(instance::DagFlowInstanceArcs, filename)
    open(filename, "w") do f
        write(f, "digraph G {\n")
        write(f, "\toverlap=scale;\n")
      	write(f, "\tsplines=true;\n")

        n,m = MyGraphModule.size(instance.g)
        for i=1:n
            println(f,"\t$i [label=\"$(instance.nodeInfoDict[i])\"]")
        end

        for i=1:length(instance.g.arcs)
            arc = instance.g.arcs[i]
            #println("(",src(e),",",dst(e),")" )
            style = ""
            if i in instance.lowArcs
                style = "style=dashed"
            end
            println(f, "\t", arc.from," -> ", arc.to, "[label = \"$i\" $style];")
        end
        write(f, "}\n")
     end
end

function writeDOTFileSol(instance::DagFlowInstanceArcs, xVal, filename)
    open(filename, "w") do f
        write(f, "digraph G {\n")
        write(f, "\toverlap=scale;\n")
      	write(f, "\tsplines=true;\n")

        n,m = MyGraphModule.size(instance.g)
        for i=1:n
            println(f,"\t$i [label=\"$(instance.nodeInfoDict[i])\"]")
        end

        for i=1:length(instance.g.arcs)
            arc = instance.g.arcs[i]
            #println("(",src(e),",",dst(e),")" )
            style = ""
            color = ""
            if i in instance.lowArcs
                style = "style=dashed"
            end
            if xVal[i] > 0.005
                color = "color=\"red\""
            end
            val = round(xVal[i],digits=2)
            println(f, "\t", arc.from," -> ", arc.to, "[label = \"$val\" $style $color];")
        end
        write(f, "}\n")
     end
end


function setupIP(instance::DagFlowInstanceArcs, file_name_instance::String)
    fileNameInstance = splitext(basename(file_name_instance))[1]
    writeDOTFile(instance, "$(fileNameInstance)_graph.dot")
    #varDict = createVarDict(g)
    g = instance.g
    n,m = MyGraphModule.size(g)
    # the nodes where flow conservation is active includes all nodes except the source and the sink:
    flowConsNodes = setdiff(1:n, [instance.sourceNode instance.sinkNode])

    IP = Model(solver=CplexSolver())

    #gInOut = createMyGraph(graph)
    @variable(IP, x[allArcs(g)] >= 0, Int)

    for i in 1:length(allArcs(g))
        arc = g.arcs[i]
        if i in instance.lowArcs
            setupperbound(x[arc], 2.0)
        else
            setupperbound(x[arc],1.0)
        end
    end



    @objective(IP, Min, sum(a.cost*x[a] for a in allArcs(g)))

    @constraint(IP, [subset in instance.arcPartition], sum(x[getArc(g,a)] for a in subset ) == 1)
    @constraint(IP, [node in flowConsNodes], sum(x[a] for a in outArcs(g,Int32(node))) == sum(x[a] for a in inArcs(g,Int32(node))))
    @constraint(IP, [outFlowLimit in instance.outFlowLimits], sum(x[a] for a in outArcs(g,outFlowLimit.node)) <= outFlowLimit.limit)
    writeLP(IP, "$(fileNameInstance)_formulation.lp", genericnames=false)

    JuMP.build(IP)
    MPB=MathProgBase
    iModel = internalmodel(IP)
    zhOrigSubSys = ZeroHalfChvatGom.findZeroHalfOrigSubSystem(MPB.getconstrmatrix(iModel), MPB.getconstrLB(iModel), MPB.getconstrUB(iModel), MPB.getvartype(iModel), false, 1)
    solve(IP, relaxation=true)
    println("objective value: ", getobjectivevalue(IP))

    xValOrig = getvalue(x)
    xVal = Vector{Float64}(undef,length(xValOrig))
    for i=1:length(xValOrig)
        xVal[i] =  xValOrig[getArc(g,i)]
        if xVal[i] > 0.001
            println("x[$i] = $(xVal[i])")
        end
    end

    strength=3
    maxIter = 40
    printCuts = true
    zhsPrintLevel = 0

    done = false
    iter = 1
    while !done && iter <= maxIter
        file_name_sol = "$(fileNameInstance)_sol_$(lpad(iter,2,"0")).dot"
        writeDOTFileSol(instance, xVal, file_name_sol)
        file_name_cuts = "$(fileNameInstance)_cuts_$(lpad(iter,2,"0")).txt"
        zhsParams = ZeroHalfChvatGom.ZeroHalfSepParams(strength,zhsPrintLevel)
        cuts = ZeroHalfChvatGom.ZeroHalfSepTranslateVars(zhOrigSubSys.ASub, zhOrigSubSys.bSub, xVal, 0.0001, zhsParams, zhOrigSubSys.varsSorted)
        f = open(file_name_cuts, "w+")
        for i = 1:length(cuts)
            #println("*** Cut $i***")
            cut = cuts[i]
            vars = cut[1]
            coeffs = cut[2]
            rhs = cut[3]
            if printCuts
                for j=1:length(vars)
                    preString = ""
                    if j > 1
                        preString = " +"
                    end
                    if coeffs[j] == 1.0
                        print(f,"$preString x$(lpad(vars[j],3,"0"))")
                    elseif coeffs[j] == -1.0
                        print(f, " - x$(lpad(vars[j],3,"0"))")
                    else
                        print(f, "$preString $(coeffs[j])*x$(lpad(vars[j],3,"0"))")
                    end
                end
                println(f, " <= $rhs")
            end
            @constraint(IP, sum(coeffs[j] * x[getArc(g,vars[j])] for j=1:length(vars)) <= rhs )
        end

        close(f)
        if length(cuts) > 0
            solve(IP, relaxation=true)

            println("Objective value (iter $iter): ", getobjectivevalue(IP))
            #println("x = ", getvalue(x))
            xValOrig = getvalue(x)
            for i=1:length(xValOrig)
                xVal[i] = xValOrig[getArc(g,i)]
            end
            iter += 1
        else
            done = true
        end
    end
end


end
