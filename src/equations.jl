function pde_odes!(dudt, u, prob, t)
    mesh = prob.geometry
    D = prob.diffusion_function
    Dp = prob.diffusion_parameters
    R = prob.reaction_function
    Rp = prob.reaction_parameters
    V = mesh.volumes
    h = mesh.spacings
    V₁ = V[begin]
    D₁ = D(u[begin], Dp)
    D₂ = D(u[begin+1], Dp)
    D̄₁₂ = (D₁ + D₂) / 2
    h₁ = h[1]
    R₁ = R(u[begin], Rp)
    dudt[begin] = inv(V₁) * (D̄₁₂ * ((u[begin+1] - u[begin]) / h₁)) + R₁
    @inbounds for i in (firstindex(dudt)+1):(lastindex(dudt)-1)
        Vᵢ = V[i]
        Dᵢ₋₁ = D(u[i-1], Dp)
        Dᵢ = D(u[i], Dp)
        Dᵢ₊₁ = D(u[i+1], Dp)
        D̄ᵢᵢ₊₁ = (Dᵢ + Dᵢ₊₁) / 2
        D̄ᵢ₋₁ᵢ = (Dᵢ₋₁ + Dᵢ) / 2
        hᵢ = h[i]
        hᵢ₋₁ = h[i-1]
        Rᵢ = R(u[i], Rp)
        dudt[i] = inv(Vᵢ) * (D̄ᵢᵢ₊₁ * ((u[i+1] - u[i]) / hᵢ) - D̄ᵢ₋₁ᵢ * ((u[i] - u[i-1]) / hᵢ₋₁)) + Rᵢ
    end
    Vₙ = V[end]
    Dₙ₋₁ = D(u[end-1], Dp)
    Dₙ = D(u[end], Dp)
    D̄ₙ₋₁ₙ = (Dₙ₋₁ + Dₙ) / 2
    hₙ₋₁ = h[end]
    Rₙ = R(u[end], Rp)
    dudt[end] = -inv(Vₙ) * (D̄ₙ₋₁ₙ * ((u[end] - u[end-1]) / hₙ₋₁)) + Rₙ
end

jacobian_sparsity(prob::FVMProblem) = jacobian_sparsity(prob.geometry.mesh_points)
function jacobian_sparsity(pts)
    num_nnz = 3(length(pts) - 2) + 4 # 3 neighbours for each interior node, 2 from each boundary node 
    I = zeros(Int64, num_nnz)
    J = zeros(Int64, num_nnz)
    V = ones(num_nnz)
    # The left boundary node 
    I[1] = firstindex(pts)
    J[1] = firstindex(pts)
    I[2] = firstindex(pts)
    J[2] = firstindex(pts) + 1
    # The interior nodes 
    ctr = 3
    for i in (firstindex(pts)+1):(lastindex(pts)-1)
        I[ctr] = i
        J[ctr] = i
        ctr += 1
        I[ctr] = i
        J[ctr] = i + 1
        ctr += 1
        I[ctr] = i
        J[ctr] = i - 1
        ctr += 1
    end
    # The right boundary node 
    I[ctr] = lastindex(pts)
    J[ctr] = lastindex(pts)
    ctr += 1
    I[ctr] = lastindex(pts)
    J[ctr] = lastindex(pts) - 1
    return sparse(I, J, V)
end