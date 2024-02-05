module CVRPInstanceModule

using Distributions
using Random
using LinearAlgebra

export CVRPInstance, plotInstance, randInstance, readInstance

mutable struct CVRPInstance
  # n: number of customers
  # node 1 to n are customers. node n+1 is the depot
  #nn (=n+1): number of nodes
  n::Int32
  nn::Int32
  Q::Int32
  coords::Array{Int64,2}
  demand::Array{Int32,1}
  # cost matrix
  c::Array{Float64,2}
  name::String
end
# default constructor
CVRPInstance() = CVRPInstance(1,2,100,[ 1 1 ], [1], [0 0; 0 0], "unknown" )

# n is the number of customers. Number of nodes is n+1 (called nn in CVRPInstance)
function randInstance(instance::CVRPInstance, n, Q, seed)
  instance.n = n
  nn =n+1
  instance.nn = nn
  instance.Q = Q
  instance.coords = zeros(Int64, nn,2)
  instance.demand = zeros(Int32, nn)
  Random.seed!(seed)
  # Generate random coordinates in the box (0,0) to (1000,1000)
  rand!(instance.coords, collect(0:1000))
  rand!(instance.demand, collect(1:30))
  # set demand of the depot to zero
  instance.demand[nn] = 0;

  # this is copy-paste from TSP example by ???
  instance.c = zeros(nn, nn)
  for i = 1:nn
      for j = i:nn
          d = norm(instance.coords[i,1:2] - instance.coords[j,1:2])
          instance.c[i,j] = d
          instance.c[j,i] = d
      end
  end
end

# n is the number of customers. Number of nodes is n+1 (called nn in CVRPInstance)
# Q: vehicle capacity
# maxDemand: customer demand is selected as a random integer in the interval [1,maxDemand]
# boxSize: customers are generated within a [0,boxSize] x [0,boxSize] box

# nCustClustered: How many customers are placed in clusters
# nClusters: How many clusters are there?
# For each cluster we generate a random  center (x',y') in the box. We also draw
# number "stdDev" as a rational number in the interval [minStdDevClust, maxStdDevClust]
# that indicates how "spread out the cluster is"
# coordinates of clusters are generated by drawing x from a truncated normal
# distribution with mean x' and std. deviation <stdDev>. The distribution is truncated to the
# interval [0,boxSize]. The same happens for the y coordinate.

# nCustConCircles: How many customers are drawn from concentric circles around the depot
# We generate <nConCircles> concentric circles. The circles are centered in midle of the
# coordinate box (independent of the depot location),
# the i'th (i=1..nConCircles) circle has radius (i/nConCircles)*(boxSize/2)
# implying that the larges has radius boxSize/2 and the smallest has radius
# (1/nConCircles)*(boxSize/2)

# depotMode: 1 = center of box, 2=one of the corners, 3=random position in the box

# @@@ to do: make demand generation "more interesting"
function randInstanceAdvanced(instance::CVRPInstance, n, Q, maxDemand, boxSize,
      nCustClustered, nClusters, minStdDevClust, maxStdDevClust,
      nCustConCircles, nConCircles,
      depotMode, seed)
  # some initial error checking:
  if nCustClustered+nCustConCircles > n
    throw("randInstanceAdvanced(...): nCustClustered + nCustConCircles > n")
  end
  if maxDemand > Q
    throw("randInstanceAdvanced(...): maxDemand > Q")
  end

  srand(seed)

  instance.n = n
  nn = n+1
  instance.nn = nn
  instance.Q = Q
  #instance.coords = zeros(Int64, nn,2)
  # generate demands:
  instance.demand = rand(1:maxDemand, nn)
  # set demand of the depot to zero
  instance.demand[nn] = 0;

  coordsUniform = []
  coordsClusters = []
  coordsCircles = []
  # Generate random coordinates in the box
  if nCustClustered+nCustConCircles < n
    nCustUniform = n - (nCustClustered + nCustConCircles)
    coordsUniform = rand(0:boxSize, nCustUniform, 2)
  end

  # Generate coordinates in clusters
  if nCustClustered > 0
    coordsClusters = zeros(Int64, nCustClustered,2)
    clusterCenters = rand(0:boxSize, nClusters, 2)
    clusterStdDev = rand(Uniform(minStdDevClust,maxStdDevClust), nClusters)

    distribX = []
    distribY = []
    for i=1:nClusters
      push!(distribX, TruncatedNormal(clusterCenters[i,1], clusterStdDev[i], 0, boxSize))
      push!(distribY, TruncatedNormal(clusterCenters[i,2], clusterStdDev[i], 0, boxSize))
    end
    for i=1:nCustClustered
      clusterId=rand(1:nClusters)
      x = round(rand(distribX[clusterId]))
      y = round(rand(distribY[clusterId]))
      coordsClusters[i,1] = x
      coordsClusters[i,2] = y
    end
  end

  # Generate coordinates in concentric circles
  if nCustConCircles > 0
    coordsCircles = zeros(Int64, nCustConCircles,2)
    radiuses = zeros(Float64, nConCircles)
    for i=1:nConCircles
      radiuses[i] = (i/nConCircles)*(boxSize/2)
    end
    for i=1:nCustConCircles
      # choose a random circle
      circleId = rand(1:nConCircles)
      # choose a random angle
      angle = rand(Uniform(0,2*pi))
      coordsCircles[i,1] = round(boxSize/2+radiuses[circleId]*cos(angle))
      coordsCircles[i,1] = round(boxSize/2+radiuses[circleId]*sin(angle))
    end
  end

  # depotMode: 1 = center of box, 2=one of the corners, 3=random position in the box
  depotCoord = zeros(Int64, 1,2)
  if depotMode == 1
    depotCoord = [boxSize/2 boxSize/2]
  elseif depotMode == 2
    depotCoord = [rand(0:1)*boxSize rand(0:1)*boxSize]
  else
    depotCoord = rand(0:boxSize, 1, 2)
  end
  instance.coords = vcat(coordsUniform, coordsClusters, coordsCircles, depotCoord)

  # this is copy-paste from TSP example by ???
  instance.c = zeros(nn, nn)
  for i = 1:nn
      for j = i:nn
          d = norm(instance.coords[i,1:2] - instance.coords[j,1:2])
          instance.c[i,j] = d
          instance.c[j,i] = d
      end
  end
end


function findKey(key::String, splitFile::Vector{SubString{String}})
    lcaseKey = lowercase(key)
    for i=1:length(splitFile)
      if lowercase(splitFile[i]) == lcaseKey
        return i
      end
    end
    return -1
end

function readInstance(fileName::String, roundDist::Bool)
    instance = CVRPInstance()
    open(fileName) do f
        #allFile = readstring(f)
        allFile = read(f, String)
        splitFile = split(allFile)

        instance.name = basename(fileName)
        # read DIMENSION
        i = findKey("DIMENSION", splitFile)
        if (i < 1)
            throw("readInstance(...): DIMENSION not found")
        end
        dim = parse(Int64, splitFile[i+2])
        instance.n = dim-1
        instance.nn = dim
        instance.coords = zeros(Int64, dim,2)
        instance.demand = zeros(Int64, dim)
        instance.c = zeros(Float64, dim, dim)
        # read VEHICLES
        # We ignore the number of vehicles specified in the file so this part is commented out
        #i = findKey("VEHICLES", splitFile)
        #=if (i >= 1)
            instance.nVeh = parse(Int32, splitFile[i+2])
        else
            instance.nVeh = 0
        end=#
        # read CAPACITY
        i = findKey("CAPACITY", splitFile)
        if (i < 1)
            throw("readInstance(...): CAPACITY not found")
        end
        instance.Q = parse(Int64, splitFile[i+2])
        # read NODE_COORD_SECTION
        i = findKey("NODE_COORD_SECTION", splitFile)
        if (i < 1)
            throw("readInstance(...): NODE_COORD_SECTION not found")
        end
        coords = zeros(Int64,dim,2)
        for j=1:dim
            index = 2+i+(j-1)*3
            x = parse(Int32, splitFile[index])
            y = parse(Int32, splitFile[index+1])
            coords[j,1] = x
            coords[j,2] = y
        end
        # move coordinates to the place we prefer (node 1:n are customers, node n+1 is the depot)
        instance.coords[dim,:] = coords[1,:]
        for j=1:dim-1
            instance.coords[j,:] = coords[j+1,:]
        end
        # compute distances
        for i = 1:dim
            for j = i:dim
                if roundDist
                    d = round(norm(instance.coords[i,1:2] - instance.coords[j,1:2]))
                else
                    d = norm(instance.coords[i,1:2] - instance.coords[j,1:2])
                end
                #println("Dist computed: $d")
                instance.c[i,j] = d
                instance.c[j,i] = d
            end
        end
        # read DEMAND_SECTION
        i = findKey("DEMAND_SECTION", splitFile)
        if (i < 1)
            throw("readInstance(...): DEMAND_SECTION not found")
        end
        demands = zeros(Int64,dim)
        for j=1:dim
            index = 2+i+(j-1)*2
            demands[j] = parse(Int32, splitFile[index])
        end
        # move demands to the place we prefer (node 1:n are customers, node n+1 is the depot)
        instance.demand[dim] = demands[1]
        for j=1:dim-1
            instance.demand[j] = demands[j+1]
        end
    end
    return instance
end

end
