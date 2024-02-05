module DAGFlowColGen

using JuMP
using CPLEX
using Gurobi
using ..MyGraphModule
using ..DagFlowInstanceModule

mutable struct RestrictedMaster
    LP::JuMP.Model
    consRef
    outflowRef
    lambda
end

function setupMaster(instance::DagFlowInstance, dummyCost::Float64)
    LP = Model(solver = CplexSolver(CPX_PARAM_SIMDISPLAY=0))
    #LP = Model(solver = GurobiSolver())

    nSets = length(instance.nodePartition)
    outFlowLimits = instance.outFlowLimits
    # Defining the variables. One per customer, corresponding to a route just containing that customer.
    @variable(LP, lambda[1:nSets] >= 0)
    # Setting the objective
    @objective(LP, Min, sum(dummyCost*lambda[i] for i=1:nSets))
    # Each partition should be "covered"
    @constraint(LP, consRef[i=1:nSets], lambda[i] == 1)
    # out flow constraints
    @constraint(LP, outflowRef[i=1:length(outFlowLimits)], 0 <= outFlowLimits[i].limit)
    RMAS = RestrictedMaster(LP, consRef, outflowRef, lambda)
    return RMAS
end

function solvePricing(instance::DagFlowInstance, dualsPartition, dualsFlowLimit, maxNPaths)
    n,m = MyGraphModule.size(instance.g)
    # pricing problem:
    nodeBonus = zeros(Float64,n)
    for i = 1:length(instance.nodePartition)
        for node = instance.nodePartition[i]
            nodeBonus[node] = dualsPartition[i]
        end
    end
    for i = 1:length(instance.outFlowLimits)
        nodeBonus[instance.outFlowLimits[i].node] += dualsFlowLimit[i]
    end
    #shortestPath(g::MyGraph, nodeOrder::Vector{Int32}, nodeBonus::Vector{Float64}, source::Int32, sink::Int32, maxNPaths::Int64, eps::Float64)
    shortestPath(instance.g, instance.nodeOrder, nodeBonus, instance.sourceNode, instance.sinkNode, maxNPaths, -0.00001)
end



function addPathsToRMAS(instance::DagFlowInstance, RMAS::RestrictedMaster, paths)

    #try to print out the nodes in the path using nodeInfoDict from instance.
    #=for i=1:length(paths)
        path = paths[i]
        println("Path $i: ")
        for node in path
            println("$node: $(instance.nodeInfoDict[node])")
        end
    end=#
    for i=1:length(paths)
        nVisits = zeros(Int64, length(instance.nodePartition))
        outFlowVisits = zeros(Int64, length(instance.outFlowLimits))
        path = paths[i]
        origCost = 0
        prevNode = -1
        for node in path
            if instance.node2Partion[node]  > 0
                nVisits[instance.node2Partion[node]] += 1
            end
            if prevNode != -1
                origCost += MyGraphModule.cost(instance.g, prevNode, node)
                #println("c($prevNode,$node) = ", MyGraphModule.cost(instance.g, prevNode, node))
            end
            if prevNode in instance.setOutFlowNodes
                for j=1:length(instance.outFlowLimits)
                    if prevNode == instance.outFlowLimits[j].node
                        outFlowVisits[j] += 1
                    end
                end
            end
            prevNode = node
        end
        #println("nVisits: $nVisits")

        touchedConstraints = ConstraintRef[]
        vals = Float64[]
        for j=1:length(nVisits)
            if nVisits[j] > 0
                push!(touchedConstraints, RMAS.consRef[j])
                push!(vals,nVisits[j])
            end
        end
        for j=1:length(instance.outFlowLimits)
            if outFlowVisits[j] > 0
                push!(touchedConstraints, RMAS.outflowRef[j])
                push!(vals,outFlowVisits[j])
            end
        end

        @variable(
            RMAS.LP,                              # Model to be modified
            xNew >= 0,                            # New variable to be added
            objective=origCost,                   # cost coefficient of new varaible in the objective
            inconstraints=touchedConstraints ,    # constraints to be modified
            coefficients=vals                     # the coefficients of the variable in those constraints
        )
        push!(RMAS.lambda, xNew) # Pushing the new variable in the array of new variables
    end
end

function colgen(instance::DagFlowInstance)
    startTime = time()
    RMAS = setupMaster(instance, 10000.0)

    iter = 0
    done = false
    timeLP = 0
    timePP = 0
    timeAddCol = 0
    colAdded = 0
    while  !done #&& iter < 100
        timeRMSStart = time()
        solve(RMAS.LP)
        timeLP += time()-timeRMSStart

        println("$iter: RMAS objective value: ", getobjectivevalue(RMAS.LP))
        #writeLP(RMAS.LP, "master.lp", genericnames=true)
        dualsPartition = getdual(RMAS.consRef)
        dualsFlowLimit = getdual(RMAS.outflowRef)
        #println("dualsPartition: $dualsPartition")
        #println("dualsFlowLimit: $dualsFlowLimit")
        timePPStart = time()
        paths, costs = solvePricing(instance, dualsPartition, dualsFlowLimit, 5)
        timePP += time()-timePPStart
        if length(paths) > 0
            timeAddColStart = time()
            addPathsToRMAS(instance, RMAS, paths)
            timeAddCol += time()-timeAddColStart
            colAdded += length(paths)
        else
            done = true
        end
        iter +=1
        println("iter = $iter, done = $done")
    end
    totTime = time()-startTime
    println("Total time: $totTime, RMAS LP time: $timeLP, PP time: $timePP, add columns time: $timeAddCol, columns added: $colAdded")

end

end
