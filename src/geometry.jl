"""
    FVMGeometry(mesh_points)
    FVMGeometry(mesh_points, spacings, volumes)

Stores the one-dimensional mesh geometry used by an [`FVMProblem`](@ref).

# Arguments

- `mesh_points`: Sorted coordinates of the finite-volume nodes.
- `spacings`: Distances between adjacent mesh points. Required only by the full constructor.
- `volumes`: Control-volume widths at the mesh points. Required only by the full
  constructor.

# Fields

- `mesh_points::T`: Sorted mesh-point coordinates.
- `spacings::T`: Distances between adjacent mesh points.
- `volumes::T`: Widths of the associated control volumes.

# Example

```julia
mesh_points = range(0.0, 1.0; length = 11)
geometry = FVMGeometry(mesh_points)
```

See also [`FVMProblem`](@ref).
"""
struct FVMGeometry{T}
    mesh_points::T
    spacings::T
    volumes::T
    function FVMGeometry(mesh_points, spacings, volumes)
        @assert issorted(mesh_points) "mesh_points is not sorted."
        @assert length(mesh_points) == length(volumes) == length(spacings) + 1 "mesh_points and volumes must have the same length, and spacings must have one less element."
        mesh_points, spacings, volumes = promote(collect(mesh_points), spacings, volumes)
        T = typeof(mesh_points)
        return new{T}(mesh_points, spacings, volumes)
    end
end
function FVMGeometry(mesh_points)
    spacings = compute_spacings(mesh_points)
    volumes = compute_volumes(mesh_points, spacings)
    return FVMGeometry(mesh_points, spacings, volumes)
end
function compute_spacings(mesh_points)
    return diff(mesh_points)
end
function compute_volumes(mesh_points, spacings)
    V = similar(mesh_points)
    V[begin] = 0.5spacings[begin]
    V[end] = 0.5spacings[end]
    for i in (firstindex(mesh_points) + 1):(lastindex(mesh_points) - 1)
        V[i] = 0.5 * (spacings[i - 1] + spacings[i])
    end
    return V
end
