module FiniteVolumeMethod1D

using SparseArrays
using SciMLBase
using CommonSolve

export FVMGeometry
export FVMProblem
export solve

include("geometry.jl")
include("problem.jl")
include("equations.jl")
include("solve.jl")

end # module 