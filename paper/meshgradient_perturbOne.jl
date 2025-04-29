
using ChainRulesCore
using FiniteDifferences
using PyCall
using Zygote
using MarineHydro

cpt = pyimport("capytaine")
resolution = (10,10)

function cptMeshPair(radius1,radius2,dx1)
    #change this to radius for pair sphere study where both are indential ; this is to only get data for getting sensitivity of one sphere with other.da11/dr2
    cptmesh1 = cpt.mesh_sphere(name="sphere", radius = radius1,center = (0,0,0),resolution = resolution).immersed_part()
    cptmesh2 = cpt.mesh_sphere(name="sphere", radius = radius2,center = (0,0,0),resolution = resolution).immersed_part()
    cptmesh = cptmesh1 + cptmesh2.translated_x(radius1+dx1) 
    return cptmesh
end


function get_center_x(J,mesh)
    centerx = mesh.centers[J,1]
    return centerx
end


function mesh_sizes(radius1,radius2,dx1) 
    cptmesh = cptMeshPair(radius1,radius2,dx1)
    return cptmesh.nb_faces, cptmesh.nb_vertices
end


h = 1e-6

∂J_r1_fd(f,r1,r2,dx1; h=1e-5) = (f(r1+h,r2, dx1) .- f(r1-h,r2, dx1)) ./ (2*h)
∂J_dx1_fd(f, r1,r2, dx1; h=1e-2) = (f(r1,r2, dx1+h) .- f(r1,r2, dx1-h)) ./ (2*h)
∂J_r2_fd(f,r1,r2,dx1; h=1e-5) = (f(r1,r2+h, dx1) .- f(r1,r2-h, dx1)) ./ (2*h)


# function define_rrule_with_finite_differences(fn, input_indices)
#     @eval function ChainRulesCore.rrule(::typeof($fn), $(Symbol.(["arg$i" for i in 1:length(input_indices)])...))
#         args = ($(Symbol.(["arg$i" for i in 1:length(input_indices)])...))
#         y = $fn($(Symbol.(["arg$i" for i in 1:length(input_indices)])...))
        
#         function f_pullback(dy)
#             df = NoTangent()  
#             gradients = Tuple(
#                 FiniteDifferences.jacobian(
#                     x -> $fn(Base.setindex(args, x, i)...), 
#                     args[i]
#                 ) for i in $input_indices
#             )
#             # Pullback: combine gradients with the adjoint `dy`
#             pullbacks = map(grad -> grad' * dy, gradients)

#             return (df, pullbacks...)
#         end

#         return y, f_pullback
#     end
# end
#mesh data from radius

function sphere_vertices(radius1,radius2,dx1) #To Do : a non differetiable argument for resolution
    cptmesh = cptMeshPair(radius1,radius2,dx1)
    return reduce(vcat,cptmesh.vertices)
end
function sphere_centers(radius1,radius2,dx1 )
    cptmesh = cptMeshPair(radius1,radius2,dx1)
    return reduce(vcat,cptmesh.faces_centers)
end

function sphere_normals(radius1,radius2,dx1)
    cptmesh = cptMeshPair(radius1,radius2,dx1)
    return reduce(vcat,cptmesh.faces_normals)
end
function sphere_radii(radius1,radius2,dx1)
    cptmesh = cptMeshPair(radius1,radius2,dx1)
    return reduce(vcat,cptmesh.faces_radiuses)
end
function sphere_areas(radius1,radius2,dx1)
    cptmesh = cptMeshPair(radius1,radius2,dx1)
    return reduce(vcat,cptmesh.faces_areas)
end

function sphere_faces(radius1,radius2,dx1)
    cptmesh = cptMeshPair(radius1,radius2,dx1)
    return reduce(vcat,cptmesh.faces)
end

# # #for one body
# define_rrule_with_finite_differences(sphere_areas, [1,2])
# define_rrule_with_finite_differences(sphere_radii, [1,2])
# define_rrule_with_finite_differences(sphere_normals, [1,2])
# define_rrule_with_finite_differences(sphere_centers, [1,2])
# define_rrule_with_finite_differences(sphere_vertices, [1,2])
# define_rrule_with_finite_differences(sphere_faces, [1,2])





"""Rules for Zygote finite differencing of the mesh"""

function ChainRulesCore.rrule(::typeof(sphere_vertices), r1,r2, dx1)
    y = sphere_vertices(r1,r2, dx1)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv_r1 = ∂J_r1_fd(sphere_vertices, r1, r2, dx1)  # Derivative w.r.t r1
        dv_r2 = ∂J_r2_fd(sphere_vertices, r1, r2, dx1)  # Derivative w.r.t r2
        dv_dx1 = ∂J_dx1_fd(sphere_vertices, r1, r2, dx1)  # Derivative w.r.t dx1
        # pullbacks
        dr1 = dv_r1' * dy
        dr2 = dv_r2' * dy
        d_dx1 = dv_dx1' * dy
        return (df, dr1, dr2, d_dx1)  
    end
    
    return y, f_pullback
end


#   function ChainRulesCore.frule(::typeof(sphere_vertices), r) #zygote by default does not use frule
#     y = NoTangent()
#     dx = ∂J_r_fd(sphere_vertices,r,dx1)' * dy

#     return y, dx
#   end




function ChainRulesCore.rrule(::typeof(sphere_centers), r1,r2,dx1)
    y = sphere_centers(r1,r2,dx1)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv_r1 = ∂J_r1_fd(sphere_centers, r1, r2, dx1)  # Derivative w.r.t r1
        dv_r2 = ∂J_r2_fd(sphere_centers, r1, r2, dx1)  # Derivative w.r.t r2
        dv_dx1 = ∂J_dx1_fd(sphere_centers, r1, r2, dx1)  # Derivative w.r.t dx1
        # pullbacks
        dr1 = dv_r1' * dy
        dr2 = dv_r2' * dy
        d_dx1 = dv_dx1' * dy
        
        return (df, dr1, dr2, d_dx1)  
    end
    return y, f_pullback
  end


function ChainRulesCore.rrule(::typeof(sphere_normals), r1,r2,dx1)
    y = sphere_normals(r1,r2,dx1)
    
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv_r1 = ∂J_r1_fd(sphere_normals, r1, r2, dx1)  # Derivative w.r.t r1
        dv_r2 = ∂J_r2_fd(sphere_normals, r1, r2, dx1)  # Derivative w.r.t r2
        dv_dx1 = ∂J_dx1_fd(sphere_normals, r1, r2, dx1)  # Derivative w.r.t dx1
        # pullbacks
        dr1 = dv_r1' * dy
        dr2 = dv_r2' * dy
        d_dx1 = dv_dx1' * dy
        return (df, dr1, dr2, d_dx1)  
    end
    return y, f_pullback
  end


function ChainRulesCore.rrule(::typeof(sphere_radii), r1,r2,dx1)
    y = sphere_radii(r1,r2,dx1)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv_r1 = ∂J_r1_fd(sphere_radii, r1, r2, dx1)  
        dv_r2 = ∂J_r2_fd(sphere_radii, r1, r2, dx1) 
        dv_dx1 = ∂J_dx1_fd(sphere_radii, r1, r2, dx1) 
        # pullbacks
        dr1 = dv_r1' * dy
        dr2 = dv_r2' * dy
        d_dx1 = dv_dx1' * dy
        return (df, dr1, dr2, d_dx1)  
    end
    return y, f_pullback
  end




function ChainRulesCore.rrule(::typeof(sphere_areas), r1,r2,dx1)
    y = sphere_areas(r1,r2,dx1)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv_r1 = ∂J_r1_fd(sphere_areas, r1, r2, dx1)  # Derivative w.r.t r1
        dv_r2 = ∂J_r2_fd(sphere_areas, r1, r2, dx1)  # Derivative w.r.t r2
        dv_dx1 = ∂J_dx1_fd(sphere_areas, r1, r2, dx1)  # Derivative w.r.t dx1
        # pullbacks
        dr1 = dv_r1' * dy
        dr2 = dv_r2' * dy
        d_dx1 = dv_dx1' * dy
        return (df, dr1, dr2, d_dx1)  
    end
    return y, f_pullback
  end



function ChainRulesCore.rrule(::typeof(sphere_faces), r1,r2,dx1)
    y = sphere_faces(r1,r2,dx1)
    function f_pullback(dy)
        df = NoTangent()  # No gradient w.r.t the function
        dv_r1 = ∂J_r1_fd(sphere_faces, r1, r2, dx1)  # Derivative w.r.t r1
        dv_r2 = ∂J_r2_fd(sphere_faces, r1, r2, dx1)  # Derivative w.r.t r2
        dv_dx1 = ∂J_dx1_fd(sphere_faces, r1, r2, dx1)  # Derivative w.r.t dx1
        # pullbacks
        dr1 = dv_r1' * dy
        dr2 = dv_r2' * dy
        d_dx1 = dv_dx1' * dy
        return (df, dr1, dr2, d_dx1)  
    end
    return y, f_pullback
  end



  function differentiableMeshPairs(r1,r2,dx1)
    nf,nv =   Zygote.@ignore mesh_sizes(r1,r2,dx1) #hack..need to know how many first to reshape
    vertices = reshape(sphere_vertices(r1,r2,dx1),(nv,3))
    faces =  reshape(sphere_faces(r1,r2,dx1),(nf,4))
    centers = reshape(sphere_centers(r1,r2,dx1),(nf,3))
    normals = reshape(sphere_normals(r1,r2,dx1),(nf,3))
    areas = sphere_areas(r1,r2,dx1)
    radii = sphere_radii(r1,r2,dx1)
    nvertices = size(vertices,1)
    nfaces = size(centers,1)
    return Mesh(vertices,faces,centers,normals,areas,radii,nvertices,nfaces) #
end

# using MarineHydro
#   #test differentiability
# function func(r1,r2,dx1)
#     mesh = differentiableMeshPairs(r1,r2,dx1)
#     S, D = assemble_matrices((Rankine(), RankineReflected(), GFWu()), mesh, 1.0)
#     return sum(imag(D)) + sum(real(S))
# end

# Zygote.jacobian(x->func(x[1],x[2],x[3]),[1.0,2.0,3.0])