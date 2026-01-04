using ..FiniteVolumeMethod1D
FVM = FiniteVolumeMethod1D

mesh_points = rand(100)
@test_throws AssertionError FVM.FVMGeometry(mesh_points)
mesh_points = sort(mesh_points)
geo = FVMGeometry(mesh_points)
@test geo.mesh_points == mesh_points
@test geo.spacings == FVM.compute_spacings(mesh_points) == diff(mesh_points)
@test geo.volumes == FVM.compute_volumes(mesh_points, geo.spacings)
@test geo.volumes ==
    0.5 *
    [geo.spacings[1]; geo.spacings[1:(end - 1)] .+ geo.spacings[2:end]; geo.spacings[end]]
