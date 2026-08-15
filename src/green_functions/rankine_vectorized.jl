# Vectorized constant-source Rankine panel from Birk (2021):
# "A Comprehensive and Practical Guide to the Hess and Smith Constant Source and Dipole Panel"
# https://doi.org/10.1080/09377255.2021.1966575
#
# Birk defines the Rankine source as φ = -1/(4πr) and evaluates potential/velocity
# in a panel-local frame (sections 2–4), then transforms vectors back with T.
# MarineHydro / Capytaine use G = 1/r; matrix assembly later multiplies by -1/(4π).
# The conversion used by VRankine is therefore
#   ∫ G dS      = -4π φ_Birk
#   ∇_{x} ∫ G dS = -4π v_Birk
# where x is the field point (element_1).

struct VRankine <: GreensFunction end
struct VRankineReflected <: GreensFunction end

const _BIRK_ZTOL = 1e-14
const _BIRK_DKTOL = 1e-14

# atan(num/den) in (-π/2, π/2) without dividing. Division by dx=0 on vertical
# edges is IEEE-finite in primal (atan(±Inf)) but Dual/Zygote NaN.
@inline function _atan_ratio(num, den)
    flipped = ifelse(den < 0, -num, num)
    return atan(flipped, abs(den))
end

@inline function _corner_svector(vertices, k)
    T = eltype(vertices)
    return SVector{3,T}(vertices[k, 1], vertices[k, 2], vertices[k, 3])
end

@inline function _corners(vertices)
    return ntuple(k -> _corner_svector(vertices, k), 4)
end

@inline function _as_center_type(center, v)
    return zero(center) .+ (v[1], v[2], v[3])
end

# Section 2.1: local panel frame from bimedians of the four corners.
function compute_local_system(corners)
    midpoints = ntuple(i -> (corners[i] + corners[mod1(i + 1, 4)]) / 2, 4)
    mgc = (midpoints[2] + midpoints[4]) / 2
    s = normalize(midpoints[2] - midpoints[4])
    t_temp = midpoints[3] - midpoints[1]
    n = normalize(cross(s, t_temp))
    t = normalize(cross(n, s))
    T = hcat(s, t, n)
    return T, mgc
end

# Section 2.2: temporary local corner coordinates with origin at the vertex centroid.
function compute_temporary_panel_corners(corners, mgc, T)
    s = T[:, 1]
    t = T[:, 2]
    Tcoord = eltype(eltype(corners))
    return ntuple(k -> SVector{3,Tcoord}(dot(corners[k] - mgc, s), dot(corners[k] - mgc, t), zero(Tcoord)), 4)
end

# Section 2.3: area and local centroid relative to the temporary origin.
function compute_area_and_centroid(temp_corners)
    T = eltype(eltype(temp_corners))
    terms = ntuple(4) do k
        xk, yk = temp_corners[k][1], temp_corners[k][2]
        xk1, yk1 = temp_corners[mod1(k + 1, 4)][1], temp_corners[mod1(k + 1, 4)][2]
        A_k = (yk1 - yk) * (xk + xk1) / 2
        Mx_k = (xk1 - xk) * (yk^2 + yk * yk1 + yk1^2) / 6
        My_k = (yk1 - yk) * (xk^2 + xk * xk1 + xk1^2) / 6
        (A_k, Mx_k, My_k)
    end
    A = sum(t -> t[1], terms)
    Mx = sum(t -> t[2], terms)
    My = sum(t -> t[3], terms)
    centroid_local = SVector{3,T}(My / A, -Mx / A, zero(T))
    return A, centroid_local
end

# Section 2.3, eq. (29): area centroid in global coordinates.
function compute_qgc(centroid_local, T, mgc)
    s = T[:, 1]
    t = T[:, 2]
    return mgc .+ centroid_local[1] .* s .+ centroid_local[2] .* t
end

# Section 2.4: local corners referred to the area centroid.
function compute_local_corners(temp_corners, centroid_local)
    return ntuple(k -> temp_corners[k] - centroid_local, 4)
end

function birk_panel_geometry(vertices)
    corners = _corners(vertices)
    T, mgc = compute_local_system(corners)
    temp_corners = compute_temporary_panel_corners(corners, mgc, T)
    area, centroid_local = compute_area_and_centroid(temp_corners)
    qgc = compute_qgc(centroid_local, T, mgc)
    local_corners = compute_local_corners(temp_corners, centroid_local)
    return T, qgc, local_corners, area
end

# Eq. (46). Uses atan (not atan2), as specified in Birk, including the paper's -1/(4π).
function velocity_potential(x, y, z, local_corners)
    contribs = ntuple(4) do k
        xk, yk = local_corners[k][1], local_corners[k][2]
        xk1, yk1 = local_corners[mod1(k + 1, 4)][1], local_corners[mod1(k + 1, 4)][2]
        dx = xk1 - xk
        dy = yk1 - yk
        dk = hypot(dx, dy)
        rk = sqrt((x - xk)^2 + (y - yk)^2 + z^2)
        rk1 = sqrt((x - xk1)^2 + (y - yk1)^2 + z^2)
        ek = (x - xk)^2 + z^2
        ek1 = (x - xk1)^2 + z^2
        hk = (x - xk) * (y - yk)
        hk1 = (x - xk1) * (y - yk1)
        log_term = log((rk + rk1 - dk) / (rk + rk1 + dk))
        log_part = ((x - xk) * dy - (y - yk) * dx) / dk * log_term
        # Rewrite (mk*ek - hk)/(z*rk) as (dy*ek - dx*hk)/(dx*z*rk) to allow vertical edges.
        atan_part = ifelse(abs(z) <= _BIRK_ZTOL, zero(z),
            z * (_atan_ratio(dy * ek - dx * hk, dx * z * rk) -
                 _atan_ratio(dy * ek1 - dx * hk1, dx * z * rk1)))
        ifelse(dk <= _BIRK_DKTOL, zero(x), -(log_part + atan_part) / 4π)
    end
    return sum(contribs)
end

# Eqs. (47)–(49), including the paper's -1/(4π). Local components are transformed
# back to global coordinates by the caller with v_g = T * v_local.
function velocity_derivatives(x, y, z, local_corners)
    contribs = ntuple(4) do k
        xk, yk = local_corners[k][1], local_corners[k][2]
        xk1, yk1 = local_corners[mod1(k + 1, 4)][1], local_corners[mod1(k + 1, 4)][2]
        dx = xk1 - xk
        dy = yk1 - yk
        dk = hypot(dx, dy)
        rk = sqrt((x - xk)^2 + (y - yk)^2 + z^2)
        rk1 = sqrt((x - xk1)^2 + (y - yk1)^2 + z^2)
        ek = (x - xk)^2 + z^2
        ek1 = (x - xk1)^2 + z^2
        hk = (x - xk) * (y - yk)
        hk1 = (x - xk1) * (y - yk1)
        log_term = log((rk + rk1 - dk) / (rk + rk1 + dk))
        vx = (dy / dk) * log_term
        vy = (-dx / dk) * log_term
        vz = ifelse(abs(z) <= _BIRK_ZTOL, zero(z),
            _atan_ratio(dy * ek - dx * hk, dx * z * rk) -
            _atan_ratio(dy * ek1 - dx * hk1, dx * z * rk1))
        ifelse(dk <= _BIRK_DKTOL, zero(SVector(x, x, x)), SVector(vx, vy, vz))
    end
    return (-1 / 4π) .* sum(contribs)
end

# Eqs. (50)–(58): Hessian of the Birk source potential in the panel-local frame.
# Only six unique entries are needed; the off-diagonals are computed independently
# so Laplace/symmetry checks can catch implementation errors. Global-frame
# Hessians are recovered with Hg = T * H * T' (eq. 37).
function velocity_hessian(x, y, z, local_corners)
    z0 = zero(x + local_corners[1][1])
    zeroH = SMatrix{3,3}(z0, z0, z0, z0, z0, z0, z0, z0, z0)
    contribs = ntuple(4) do k
        xk, yk = local_corners[k][1], local_corners[k][2]
        xk1, yk1 = local_corners[mod1(k + 1, 4)][1], local_corners[mod1(k + 1, 4)][2]
        dx = xk1 - xk
        dy = yk1 - yk
        dk = hypot(dx, dy)
        if dk <= _BIRK_DKTOL
            return zeroH
        end
        dxk = x - xk
        dyk = y - yk
        dxk1 = x - xk1
        dyk1 = y - yk1
        rk = sqrt(dxk^2 + dyk^2 + z^2)
        rk1 = sqrt(dxk1^2 + dyk1^2 + z^2)
        rsum = rk + rk1
        den = rsum^2 - dk^2
        ρk = rk * rk1 + dxk * dxk1 + dyk * dyk1 + z^2
        λk = dxk * dyk1 - dxk1 * dyk
        xterm = dxk / rk + dxk1 / rk1
        yterm = dyk / rk + dyk1 / rk1
        zterm = z / rk + z / rk1
        ρden = (rk * rk1) * ρk
        hxx = 2 * dy / den * xterm
        hyx = -2 * dx / den * xterm
        hxy = 2 * dy / den * yterm
        hyy = -2 * dx / den * yterm
        hxz = 2 * dy / den * zterm
        hyz = -2 * dx / den * zterm
        hzx = z * dy * rsum / ρden
        hzy = -z * dx * rsum / ρden
        hzz = λk * rsum / ρden
        # Column-major SMatrix: columns are (∂/∂x, ∂/∂y, ∂/∂z) of (φx, φy, φz).
        SMatrix{3,3}(hxx, hyx, hzx, hxy, hyy, hzy, hxz, hyz, hzz)
    end
    return (-1 / 4π) .* sum(contribs)
end

function _flatten_corners(local_corners)
    return SVector{12}(
        local_corners[1][1], local_corners[1][2], local_corners[1][3],
        local_corners[2][1], local_corners[2][2], local_corners[2][3],
        local_corners[3][1], local_corners[3][2], local_corners[3][3],
        local_corners[4][1], local_corners[4][2], local_corners[4][3],
    )
end

function _unflatten_corners(v)
    return (
        SVector(v[1], v[2], v[3]),
        SVector(v[4], v[5], v[6]),
        SVector(v[7], v[8], v[9]),
        SVector(v[10], v[11], v[12]),
    )
end

# Geometry (panel corners) has no closed-form Birk derivative; fill those
# cotangents with ForwardDiff through the primal. Duals take the primal path,
# so these pullbacks do not recurse into the rrules below.
function _potential_corners_adjoint(ȳ, x, y, z, local_corners)
    g = ForwardDiff.gradient(_flatten_corners(local_corners)) do v
        return velocity_potential(x, y, z, _unflatten_corners(v))
    end
    return _unflatten_corners(ȳ .* g)
end

function _velocity_corners_adjoint(ȳ, x, y, z, local_corners)
    J = ForwardDiff.jacobian(_flatten_corners(local_corners)) do v
        return velocity_derivatives(x, y, z, _unflatten_corners(v))
    end
    return _unflatten_corners(J' * ȳ)
end

function _cotangent_svector(ȳ)
    ȳu = unthunk(ȳ)
    return SVector{3}(ȳu[1], ȳu[2], ȳu[3])
end

# Zygote / ChainRules: dφ/d(x,y,z) is the analytic Birk velocity (eqs. 47–49).
function ChainRulesCore.rrule(::typeof(velocity_potential), x, y, z, local_corners)
    φ = velocity_potential(x, y, z, local_corners)
    v = velocity_derivatives(x, y, z, local_corners)
    function velocity_potential_pullback(ȳ)
        if ȳ isa AbstractZero
            return (NoTangent(), ZeroTangent(), ZeroTangent(), ZeroTangent(), ZeroTangent())
        end
        ȳu = unthunk(ȳ)
        Δcorners = _potential_corners_adjoint(ȳu, x, y, z, local_corners)
        return (NoTangent(), ȳu * v[1], ȳu * v[2], ȳu * v[3], Δcorners)
    end
    return φ, velocity_potential_pullback
end

# Zygote / ChainRules: dv/d(x,y,z) is the analytic Birk Hessian (eqs. 50–58).
function ChainRulesCore.rrule(::typeof(velocity_derivatives), x, y, z, local_corners)
    v = velocity_derivatives(x, y, z, local_corners)
    H = velocity_hessian(x, y, z, local_corners)
    function velocity_derivatives_pullback(ȳ)
        if ȳ isa AbstractZero
            return (NoTangent(), ZeroTangent(), ZeroTangent(), ZeroTangent(), ZeroTangent())
        end
        ȳv = _cotangent_svector(ȳ)
        gx = H' * ȳv
        Δcorners = _velocity_corners_adjoint(ȳv, x, y, z, local_corners)
        return (NoTangent(), gx[1], gx[2], gx[3], Δcorners)
    end
    return v, velocity_derivatives_pullback
end

function _vrankine_integral_and_gradient(element_1, element_2)
    point = center(element_1)
    source_point = center(element_2)
    source_radius = radius(element_2)
    source_area = area(element_2)
    r̄ = point - source_point
    r = norm(r̄)
    if r > 7 * source_radius
        # Same far-field kernel as Rankine: ∫ G dS ≈ A/r and ∇_{x2} ∫ G ≈ (x1-x2) A / r^3.
        return source_area / r, (r̄ / r^3) * source_area
    end

    Tmat, qgc, local_corners, _area = birk_panel_geometry(vertices(element_2))
    PT = eltype(point)
    QT = eltype(qgc)
    field_local = Tmat' * (SVector{3,PT}(point[1], point[2], point[3]) - SVector{3,QT}(qgc[1], qgc[2], qgc[3]))
    φ_birk = velocity_potential(field_local[1], field_local[2], field_local[3], local_corners)
    v_local = velocity_derivatives(field_local[1], field_local[2], field_local[3], local_corners)
    v_global = Tmat * v_local
    # Convert Birk (φ = -1/(4πr)) to MarineHydro (G = 1/r).
    φ = -4π * φ_birk
    grad_wrt_field = -4π .* v_global   # ∇_{x1} ∫ G dS
    # Rankine far/near API stores the gradient-with-respect-to-second-variable
    # as the default `integral_gradient` value (see rankine.jl).
    grad_wrt_source = -grad_wrt_field
    return φ, grad_wrt_source
end

function greens(::VRankine, element_1, element_2, wavenumber=nothing)
    return greens(Rankine(), element_1, element_2, wavenumber)
end

function gradient_greens(::VRankine, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)
    return gradient_greens(Rankine(), element_1, element_2, wavenumber; with_respect_to_first_variable)
end

function integral(::VRankine, element_1, element_2, wavenumber=nothing)
    φ, _ = _vrankine_integral_and_gradient(element_1, element_2)
    return φ
end

function integral_gradient(::VRankine, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)
    _, grad_wrt_source = _vrankine_integral_and_gradient(element_1, element_2)
    point = center(element_1)
    if with_respect_to_first_variable
        return _as_center_type(point, -grad_wrt_source)
    else
        return _as_center_type(point, grad_wrt_source)
    end
end

function both_integral_and_integral_gradient(::VRankine, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)
    φ, grad_wrt_source = _vrankine_integral_and_gradient(element_1, element_2)
    point = center(element_1)
    if with_respect_to_first_variable
        return φ, _as_center_type(point, -grad_wrt_source)
    else
        return φ, _as_center_type(point, grad_wrt_source)
    end
end

function greens(::VRankineReflected, element_1, element_2, wavenumber=nothing)
    return greens(VRankine(), free_surface_symmetry(element_1), element_2)
end

function gradient_greens(::VRankineReflected, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)
    if with_respect_to_first_variable
        return gradient_greens(VRankine(), element_1, free_surface_symmetry(element_2); with_respect_to_first_variable)
    else
        return gradient_greens(VRankine(), free_surface_symmetry(element_1), element_2; with_respect_to_first_variable)
    end
end

function integral(::VRankineReflected, element_1, element_2, wavenumber=nothing)
    return integral(VRankine(), free_surface_symmetry(element_1), element_2)
end

function integral_gradient(::VRankineReflected, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)
    ng = integral_gradient(VRankine(), free_surface_symmetry(element_1), element_2; with_respect_to_first_variable)
    if with_respect_to_first_variable
        return vertical_reflection(ng)
    else
        return ng
    end
end

function both_integral_and_integral_gradient(::VRankineReflected, element_1, element_2, wavenumber=nothing; with_respect_to_first_variable=false)
    g, ng = both_integral_and_integral_gradient(VRankine(), free_surface_symmetry(element_1), element_2; with_respect_to_first_variable)
    if with_respect_to_first_variable
        return g, vertical_reflection(ng)
    else
        return g, ng
    end
end
