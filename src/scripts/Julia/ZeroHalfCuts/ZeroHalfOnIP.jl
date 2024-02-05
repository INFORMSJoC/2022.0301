include("ZeroHalfCuts.jl")

module ZeroHalfOnIP

using JuMP
using CPLEX
using LightGraphs
using SimpleWeightedGraphs
using SparseArrays

using ..ZeroHalfChvatGom
using ..readLP

function writeDOTFile(g::SimpleWeightedDiGraph, filename)
    open(filename, "w") do f
        w = weights(g)
        write(f, "digraph G {\n")
        write(f, "\toverlap=scale;\n")
      	write(f, "\tsplines=true;\n")

        n = nv(g)
        for e in edges(g)
            #println("(",src(e),",",dst(e),")" )
            println(f, "\t",src(e)," -> ", dst(e), "[",
                    "label = \"", w[src(e),dst(e)], "\"];")
        end
        write(f, "}\n")
     end
end

function findGraph(A,consLB, consUB)
    # go through the constraints. Any ... = 0 constraint will give rise to a node in the graph
    (m,n) = size(A)
    Atranspose = sparse(A')
    nNodes = 0
    # each variable corresponds to an arc.
    # the head an tail variables define each arc
    head = zeros(Int64, n)
    tail = zeros(Int64, n)

    rows = rowvals(Atranspose)
    vals = nonzeros(Atranspose)
    for i=1:m
        if consLB[i] == consUB[i] && consLB[i] != 1
            #println("Constraint $i ... ")
            nNodes += 1
            # go through  row i of A (the same as column i in A^T )
            for j in nzrange(Atranspose, i)
                varID = rows[j]
                val = vals[j]
                if val != 1 && val != 0 && val != -1
                    throw("Unexpected value in flow assumed conservation constraint")
                else
                    if val < 0
                        if tail[varID] != 0
                            throw("constraint $i: arc $varID already has a tail")
                        end
                        tail[varID] = nNodes
                    elseif val > 0
                        if head[varID] != 0
                            throw("constraint $i: arc $varID already has a head")
                        end
                        head[varID] = nNodes
                    end
                end
            end
        end
    end
    g = SimpleWeightedDiGraph(nNodes+2)
    for i=1:n
        if tail[i] == 0
            throw("arc $i has no tail")
        elseif head[i] == 0
            throw("arc $i has no head")
        else
            add_edge!(g, tail[i], head[i], i)
        end
    end
    writeDOTFile(g, "graph.dot")
end

function test(lpFileName, maxIter=1, strength=3, printCuts=false)
    A, obj, varLB, varUB, consLB, consUB, varType = readLP.loadLP(lpFileName)
    findGraph(A,consLB, consUB)
    zhOrigSubSys = ZeroHalfChvatGom.findZeroHalfOrigSubSystem(A, consLB, consUB, varType, false, 1)
    (m,n) = size(A)
    model = Model(solver = CplexSolver(CPX_PARAM_SCRIND=0))
    @variable(model, varLB[j] <= x[j=1:n] <= varUB[j] )

    @objective(model, Min, sum( obj[j] * x[j] for j=1:n) )
    @constraint(model, consLB .<= A*x .<= consUB )
    solve(model)

    println("Objective value: ", getobjectivevalue(model))
    #println("x = ", getvalue(x))
    xVal = getvalue(x)
    for i=1:length(xVal)
        if xVal[i] > 0.001
            println("x[$i] = $(xVal[i])")
        end
    end
    #writeLP(model,"test.lp")

    done = false
    iter = 1
    while !done && iter <= maxIter
        zhsParams = ZeroHalfChvatGom.ZeroHalfSepParams(strength,1)
        cuts = ZeroHalfChvatGom.ZeroHalfSepTranslateVars(zhOrigSubSys.ASub, zhOrigSubSys.bSub, xVal, 0.0001, zhsParams, zhOrigSubSys.varsSorted)
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
                        print("$preString x$(vars[j])")
                    elseif coeffs[j] == -1.0
                        print(" - x$(vars[j])")
                    else
                        print("$preString $(coeffs[j])*x$(vars[j])")
                    end
                end
                println(" <= $rhs")
            end
            @constraint(model, sum(coeffs[j] * x[vars[j]] for j=1:length(vars)) <= rhs )
        end
        if length(cuts) > 0
            solve(model)

            println("Objective value (iter $iter): ", getobjectivevalue(model))
            #println("x = ", getvalue(x))
            xVal = getvalue(x)
            iter += 1
        else
            done = true
        end
    end
end

ZeroHalfOnIP.test("LP_and_blocks/Misc/wt008_028_2_mip.lp")
#ZeroHalfOnIP.test("LP_and_blocks/Misc/DagIntFlow_A-n32.lp")

end
