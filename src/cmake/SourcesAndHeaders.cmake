set(sources
    src/BddCoeff.cpp
    src/BranchBoundTree.cpp
    src/BranchNode.cpp
    src/Column.cpp
    src/IO.cpp
    src/Instance.cpp
    src/Interval.cpp
    src/Job.cpp
    src/LocalSearch.cpp
    src/Lowerbound.cpp
    src/Model.cpp
    src/ModelInterface.cpp
    src/NodeBdd.cpp
    src/NodeData.cpp
    src/Parms.cpp
    src/PricerSolverArcTimeDP.cpp
    src/PricerSolverBase.cpp
    src/PricerSolverBdd.cpp
    src/PricerSolverBddBackward.cpp
    src/PricerSolverBddForward.cpp
    src/PricerSolverSimpleDP.cpp
    src/PricerSolverWrappers.cpp
    src/PricingStabilization.cpp
    src/Problem.cpp
    src/SeperationSolver.cpp
    src/Solution.cpp
    src/StabilizationWrappers.cpp
    src/Statistics.cpp
    src/VariableKeyBase.cpp
    # src/ZeroHalfCuts.cpp
    # src/ZeroHalfSystem.cpp
)

set(exe_sources "src/main.cpp" ${sources})

set(headers ${CMAKE_CURRENT_SOURCE_DIR}/include)
