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

@static if VERSION < v"1.7"
    struct Returns{V} <: Function
        value::V
        Returns{V}(value) where {V} = new{V}(value)
        Returns(value) = new{Base._stable_typeof(value)}(value)
    end
    (obj::Returns)(@nospecialize(args...); @nospecialize(kw...)) = obj.value
end

end # module 