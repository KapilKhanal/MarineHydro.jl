
using Plots, ColorTypes
using CSV
using LaTeXStrings
using DataFrames, Statistics

orange = RGB(230/255,159/255,0/255)  
vermillion = RGB(213/255, 94/255, 0/255) 
bluishgreen = RGB(0/255, 158/255, 115/255)  
default(size = (800, 600))
default(
    guidefont = font(15),   # Axis labels
    tickfont = font(12),    # Tick labels
    legendfont = font(10)   # Legend text
)

data = CSV.File("/home/cornell/ForkMarineHydro/MarineHydro.jl/paper/Plots/12_heuristics_dx_DELhommeau_singleperturb.csv") |> DataFrame
omega = 1.03
k = omega^2/9.8

period = 2*pi*9.8/(k^2)


plot( data.dx, data.A12_grad_r, label=L"\frac{\partial A_{12}(\omega = 1.03\,\mathrm{rad/s})}{\partial r_2}", marker=:circle, lw=2, color=bluishgreen,legend=:right)
plot!(data.dx, data.B12_grad_r, label=L"\frac{\partial B_{12}(\omega = 1.03\,\mathrm{rad/s})}{\partial r_2}", marker=:square, lw=2, color=vermillion)
vline!( [10], label="PWA heuristics\n (Singh and Babarit (2013))", lw=2, linestyle=:dash, color=orange)

# Labels and Title
xlabel!("kx",fontsize=18)
ylabel!("Sensitivity values",fontsize=18)

# Save the plot
savefig("/home/cornell/ForkMarineHydro/MarineHydro.jl/paper/Plots/12_heuristics_dx_r2.pdf")


# # #switch for damping and added mass here and corresponding labels , columns below
# data = CSV.File("/home/cornell/BEMJulia/MarineHydro.jl/paper/Plots/added_mass_data_dimensionless.csv") |> DataFrame

# # Unique values for dimensionless parameters
# dx_r_ratios = unique(data.dx_r_ratio)
# kr_values = unique(data.kr)

# dx_values = unique(data.dx_r_ratio)
# kr_values =  unique(data.kr)

# grad_r_matrix = reshape(data.grad_r, length(dx_values), length(kr_values))

# # Normalize function
# function normalize(matrix)
#     min_val = minimum(matrix)
#     max_val = maximum(matrix)
#     println("Maximum: $max_val")
#     println("Minimum: $min_val")
#     return (matrix .- min_val) ./ (max_val - min_val)
# end

# # Normalize gradient matrices
# grad_r_matrix = normalize(grad_r_matrix)


# # Plot heatmap using dimensionless parameters
# p1 = heatmap(
#     dx_r_ratios, kr_values, grad_r_matrix,
#     ylabel="Kr",
#     xlabel="x/r",
#     color=:magma,
#     xguidefontsize=24,
#     yguidefontsize=24,
#     titlefontsize=16,
#     tickfontsize=8,
#     colorbar_title=L"\frac{\partial A_{11}(\omega = 1.03\,\mathrm{rad/s})}{\partial r_2}", #switch
#     colorbar_titlefontsize=15,
#     colorbar_titleorientation=:vertical
# )

# savefig("//home/cornell/ForkMarineHydro/MarineHydro.jl/paper/Plots/added_mass_dimensionless_grad_dr.pdf")

