using Plots
using CSV
using DataFrames
using MarineHydro
include("/home/cornell/ForkMarineHydro/MarineHydro.jl/paper/meshgradient_perturbOne.jl")
 #takes a while depending on this - some faces_max_radius gives weird answer
#check meshes.jl to change faces_max_radius.
radius_range = [1,2,3,4,5]
heave = [0,0,1] #heave
surge = [1,0,0]

plot(size=(800, 600))



#impact of one body on the other added mass coeffficients
function added_mass_off_diagonal(radius1,radius2,omega ,dx1)  
    mesh = differentiableMeshPairs(radius1,radius2, dx1)  
    total_nfaces = mesh.nfaces 
    element_is_in_sphere_1(j) = get_center_x(j, mesh) <= radius1
    element_is_in_sphere_2(j) = get_center_x(j, mesh) > radius1
    sphere_1_heave_normal = [element_is_in_sphere_1(j) ? mesh.normals[j,:]' * [0.0, 0.0, 1.0] : 0.0 for j in 1:total_nfaces]
    sphere_2_heave_normal = [element_is_in_sphere_2(j) ? mesh.normals[j,:]' * [0.0, 0.0, 1.0] : 0.0 for j in 1:total_nfaces]
    k = omega^2 / 9.81  # Wave number
    S, D = assemble_matrices((Rankine(), RankineReflected(), ExactGuevelDelhommeau()), mesh, k)
    potential = MarineHydro.solve(D, S, -1im * omega * sphere_1_heave_normal)
    pressure = 1im * 1000 * omega * potential
    # force = -sum(pressure .* sphere_2_heave_normal .* mesh.areas)
    # A12 = real(force) / omega^2
    force = -sum(pressure .* sphere_1_heave_normal .* mesh.areas)
    A11 = real(force) / omega^2
    return A11
end

function damping_off_diagonal(radius1,radius2,omega ,dx1)  
    mesh = differentiableMeshPairs(radius1,radius2, dx1)  
    total_nfaces = mesh.nfaces 
    element_is_in_sphere_1(j) = get_center_x(j, mesh) <= radius1
    element_is_in_sphere_2(j) = get_center_x(j, mesh) > radius1
    sphere_1_heave_normal = [element_is_in_sphere_1(j) ? mesh.normals[j,:]' * [0.0, 0.0, 1.0] : 0.0 for j in 1:total_nfaces]
    sphere_2_heave_normal = [element_is_in_sphere_2(j) ? mesh.normals[j,:]' * [0.0, 0.0, 1.0] : 0.0 for j in 1:total_nfaces]
    k = omega^2 / 9.81  # Wave number
    S, D = assemble_matrices((Rankine(), RankineReflected(), ExactGuevelDelhommeau()), mesh, k) # Assemble matrices tuple error -- use default
    potential = MarineHydro.solve(D, S, -1im * omega * sphere_1_heave_normal)
    pressure = 1im * 1000 * omega * potential
    force = -sum(pressure .* sphere_1_heave_normal .* mesh.areas)
    B11 = imag(force) / omega
    # force = -sum(pressure .* sphere_2_heave_normal .* mesh.areas)
    # B12 = imag(force) / omega
    return B11
end

# test zygote computed with finite difference
function finite_diff_grad_r2(f, r1, r2, ω, dx1; ϵ=1e-6)
    r2_plus = r2 + ϵ
    r2_minus = r2 - ϵ
    (f(r1, r2_plus, ω, dx1) - f(r1, r2_minus, ω, dx1)) / (2ϵ)
end

# r1, r2, ω = 1.0, 1.0, 1.03
# dx1 = 2.0
# fd_grad = finite_diff_grad_r2(added_mass_off_diagonal, r1, r2, ω, dx1)  #w.r.t r2
# zygote_grad = Zygote.gradient(r2 -> added_mass_off_diagonal(r1,r2,ω,dx1),  r2)

# println("Finite difference gradient: ", fd_grad)
# println("Zygote gradient: ", zygote_grad)

using MarineHydro
# Set parameters --change ii,ij depending on which entry of A,B
g = 9.8 
heave = [0, 0, 1]  # Heave
dx_r_ratios = collect(range(2.0, stop=200.5, step=4))
kr_values = collect(range(0.5, stop=12.0, step=1.0))  # k*r dimensionless parameter
r1 = 1.0   
omega = 1.03
r2 = 1.0

# # #check with heuristics that the  - plot sensitivity of added mass and damping with dx and see if they go to zero or low.
# data = DataFrame(A11_grad_r=Float64[], B11_grad_r =Float64[], dx = Float64[], A11=Float64[], B11=Float64[])
# for dx in dx_r_ratios
#     @show dx
#     A11_grad_r, = Zygote.gradient(r2 -> added_mass_off_diagonal(r1,r2,omega ,dx), r2) #at fix r = 1.0
#     B11_grad_r, =  Zygote.gradient(r2 -> damping_off_diagonal(r1,r2,omega ,dx), r2)
#     @show A11_grad_r
#     @show B11_grad_r
#     A11 = added_mass_off_diagonal(r1,r2,omega ,dx)
#     B11 = damping_off_diagonal(r1,r2,omega ,dx)
#     push!(data, (A11_grad_r, B11_grad_r, dx, A11, B11))
# end
# print(data)
# CSV.write("/home/cornell/ForkMarineHydro/MarineHydro.jl/paper/Plots/more_11_heuristics_dx_DELhommeau_singleperturb.csv", data)


dx_r_ratios = collect(range(1.5, stop=20.5, step=0.5))
kr_values = collect(range(0.5, stop=12.0, step=0.1))  # k*r dimensionless parameter
data = DataFrame(dx_r_ratio=Float64[], kr=Float64[], grad_r = Float64[])

for dx_r in dx_r_ratios
    for kr in kr_values
        dx1 = dx_r * r1
        @show dx1
        omega = sqrt(kr * g / r1)
        @show omega
        grad_r, = Zygote.gradient(r2 -> added_mass_off_diagonal(r1,r2,omega ,dx1), r2)
        @show grad_r
        push!(data, (dx_r, kr, grad_r))
    end
end

CSV.write("/home/cornell/ForkMarineHydro/MarineHydro.jl/paper/Plots/single_perturb_added_mass_data_dimensionless.csv", data)

println("damping data")
data = DataFrame(dx_r_ratio=Float64[], kr=Float64[], grad_r = Float64[])

for dx_r in dx_r_ratios
    for kr in kr_values
        dx1 = dx_r * r1
        omega = sqrt(kr * g / r1)
        grad_r, =  Zygote.gradient(r2 -> damping_off_diagonal(r1,r2,omega ,dx1), r2)
        @show grad_r
        push!(data, (dx_r, kr, grad_r))
    end
end

CSV.write("/home/cornell/ForkMarineHydro/MarineHydro.jl/paper/Plots/single_perturb_damping_data_dimensionless.csv", data)
