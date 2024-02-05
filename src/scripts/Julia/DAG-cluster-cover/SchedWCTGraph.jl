module SchedWCTGraph

mutable struct SchedWCTinstance
    m::Int64 # number of machines
    weight::Vector{Int64}
    procTime::Vector{Int64}
    releaseTime::Vector{Int64}
    deadline::Vector{Int64}
end

SchedWCTinstance() = SchedWCTinstance(0,[],[],[],[])

function readSchedWCTInstance(fileName::String, m::Int64)
    instance = SchedWCTinstance()
    open(fileName) do f
        allFile = read(f, String)
        splitFile = split(allFile)
        #println(splitFile)
        n = parse(Int64,splitFile[2])
        weight = zeros(Int64, n)
        procTime = zeros(Int64, n)
        releaseTime = zeros(Int64, n)
        deadline = typemax(Int64)*ones(Int64,n)
        for i=1:n
            #println("i=$i")
            weight[i] = parse(Int64,splitFile[1+3*i])
            procTime[i] = parse(Int64,splitFile[2+3*i])
        end
        instance = SchedWCTinstance(m,weight, procTime, releaseTime, deadline)
    end
    return instance
end


function constructSchedWCTGraph(instance::SchedWCTinstance)
    n = length(instance.weight)
    m = instance.m
    itemOrder = sort(collect(1:n), by = i -> (instance.weight[i] / instance.procTime[i]), rev=true)
    #itemOrder = sort(collect(1:n))
    println("itemOrder: ", itemOrder)
    for i in itemOrder
        println("$i: weight: $(instance.weight[i]), processing time: $(instance.procTime[i]), ratio: $(instance.weight[i] / instance.procTime[i])")
    end
    # compute latest completion time based on property #3 from A Branch-and-Price Algorithm for Parallel Machine
    # Scheduling Using ZDDs and Generic Branching (Daniel Kowalczyk, Roel Leus)
    HMax = sum(instance.procTime) / m + (m-1)/m*maximum(instance.procTime)
    println("HMax = $HMax")
end

end

#SchedWCTGraph.readSchedWCTInstance("Data/scheduling_WCT/instances_wct/20_3/instance1_20_3_0.txt")
#SchedWCTGraph.constructSchedWCTGraph(x)
