struct DOF
    name::Symbol
    direction::Vector{Float64} 
end

"""
Individual body class with its name, mesh, and degrees of freedom.

# Fields
- `mesh::Mesh`: The mesh geometry
- `name::String`: Body name
- `dofs::Matrix{Float64}`: Degrees of freedom matrix (ndofs × 3) where each row is a DOF vector
"""
struct FloatingBody
    mesh::Mesh
    name::String
    dofs::Vector{DOF}
    
    function FloatingBody(mesh::Mesh, name::String, dofs::Vector{DOF})
        return new(mesh, name, dofs)
    end
    #if no dofs are provided, create an empty vector of DOFs
    function FloatingBody(mesh::Mesh)
        println("no dofs provided for body")
        return new(mesh, "unnamed_body", DOF[])
    end
end

"""
A type representing a collection of floating bodies for multi-body multi-dof hydrodynamics computations.

# Fields
- `bodies::Vector{FloatingBody}`: The floating bodies
- `mesh::Mesh`: The combined mesh of all bodies
- `panel_to_body::Vector{Int}`: Mapping from panel index to body index
"""
struct MultiBody
    mesh::Mesh
    body_panels::Vector{UnitRange{Int}}
    dofs::Vector{Vector{DOF}} #each body has a vector of DOFs
end



"""
    MultiBody(bodies::Vector{FloatingBody})
"""
function MultiBody(bodies::Vector{FloatingBody})
    meshes = [body.mesh for body in bodies]
    panelCounts = [m.nfaces for m in meshes]
    # Concatenate all mesh fields?
    vertices = vcat([m.vertices for m in meshes]...)
    faces = vcat([m.faces for m in meshes]...)
    centers = vcat([m.centers for m in meshes]...)
    normals = vcat([m.normals for m in meshes]...)
    areas = vcat([m.areas for m in meshes]...)
    radii = vcat([m.radii for m in meshes]...)
    nvertices = sum(m.nvertices for m in meshes)
    nfaces = sum(panelCounts)
    dofs = [body.dofs for body in bodies]
    # combined mesh
    starts = cumsum([1; panelCounts[1:end-1]])
    stops = cumsum(panelCounts)
    body_panels = [s:e for (s, e) in zip(starts, stops)]
    println("body_panels: ", body_panels)
    mesh = Mesh(vertices, faces, centers, normals, areas, radii, nvertices, nfaces)
    return MultiBody(mesh, body_panels, dofs)
end

"""
    compute_multibody_radiation(mb::MultiBody, ω::Real)

Compute the full added mass and damping matrices for a MultiBody system at frequency ω, for arbitrary DOFs per body.

# Arguments
- `mb`: MultiBody object
- `ω`: Angular frequency
- `dofs`: Vector of vectors, each containing the DOFs for a body (e.g., `[[1,0,0], [0,0,1]]`)

# Returns
- `A`: Added mass matrix (nbodies, ndof_k, nbodies, ndof_l)
- `B`: Damping matrix (nbodies, ndof_k, nbodies, ndof_l)
"""
function compute_multibody_radiation(mb::MultiBody, ω::Real)
    nbodies = length(mb.body_panels)
    # each body can have multiple DOFs
    dofs = [mb.dofs[i] for i in 1:nbodies]  
    ndofs = [length(dofs[i]) for i in 1:nbodies] 
    maxdof = maximum(ndofs) # maximum number of DOFs for any body that interacts with other bodies
    
    A = zeros(Float64, nbodies, maxdof, nbodies, maxdof)
    B = zeros(Float64, nbodies, maxdof, nbodies, maxdof)
    forces = zeros(ComplexF64, nbodies, maxdof, nbodies, maxdof)
    
    mesh = mb.mesh
    ρ = 1000.0
    
    
    # Precompute generalized normal vectors for all DOFs of all bodies
    # G[k, i, p] = n_i^k for panel p of body k and dof i
    G = zeros(Float64, nbodies, maxdof, mesh.nfaces) 
    
    for k in 1:nbodies
        panels_k = mb.body_panels[k]
        for i in 1:ndofs[k]
            for p in panels_k
                G[k, i, p] = -dot(dofs[k][i].direction, mesh.normals[p, :])
            end
        end
    end
    
    # For each dof of each body (body l, DOF j)
    for l in 1:nbodies
        for j in 1:ndofs[l] #ndof has counts for the dofs for each body 
            # Build radiation BC:  body l with DOF j
            bc = zeros(ComplexF64, mesh.nfaces)
            panels_l = mb.body_panels[l]
            println("panels_l: ", panels_l)
            bc[panels_l] = radiation_bc(mb.mesh, dofs[l][j].direction, ω)[panels_l] #messy needs individual mesh 
            
            # Solve for potential Φ_j^l across all bodies (mesh faces)
            k = ω^2 / 9.8
            S, D = assemble_matrices((Rankine(), RankineReflected(), GFWu()), mesh, k)
            potential = MarineHydro.solve(D, S, bc)
            
            for k in 1:nbodies
                panels_k =  mb.body_panels[k]
                potential_k = potential[panels_k]
                areas_k = mesh.areas[panels_k]
                
                for i in 1:ndofs[k]
                    # generalized normal for this DOF and body
                    # G[k, i, 1:length(panels_k)] contains the normal amplitudes
                    normal_dof_amp = G[k, i, 1:length(panels_k)]
                    
                    # F_ij^(k)(l) = -ρ ∬_S_k Φ_j^l n_i^k dS
                    #  -ρ * (normal_dof_amp .* areas_k)' * potential_k
                    f = -ρ * dot(normal_dof_amp .* areas_k, potential_k)
                    
                    forces[k, i, l, j] = f
                    A[k, i, l, j] = real(f) / ω^2
                    B[k, i, l, j] = imag(f) / ω
                end
            end
        end
    end
    return A, B
end

# """
#     compute_multibody_diffraction(mb::MultiBody, ω::Real, dofs::Vector{Vector{<:Real}})

# Compute the Froude-Krylov and diffraction forces for a MultiBody system at frequency ω, for arbitrary DOFs per body.

# # Arguments
# - `mb`: MultiBody object
# - `ω`: Angular frequency
# - `dofs`: Vector of vectors, each containing the DOFs for a body (e.g., `[[1,0,0], [0,0,1]]`)

# # Returns
# - `fk_forces`: Froude-Krylov forces (nbodies, maxdof)
# - `diff_forces`: Diffraction forces (nbodies, maxdof)
# """
# function compute_multibody_diffraction(mb::MultiBody, ω::Real, dofs::Vector{Vector{<:Real}})
#     nbodies = length(mb.bodies)
#     ndofs = [length(d) for d in dofs]
#     maxdof = maximum(ndofs)
#     fk_forces = zeros(Float64, nbodies, maxdof)
#     diff_forces = zeros(Float64, nbodies, maxdof)
#     for b in 1:nbodies
#         for j in 1:ndofs[b]
#             fk, diff = solve_diffraction_problem(mb.bodies[b], ω, dofs[b][j])
#             fk_forces[b, j] = fk
#             diff_forces[b, j] = diff
#         end
#     end
#     return fk_forces, diff_forces
# end


