function As(x) 
	if 0 <= x <= .5
		return 1 + 2.2(x - 0.5)^2
	elseif .5 <= x <= 1
		return 1 + .22(x - 0.5)^2
	end
end

function AB2(z₀, x, t, t_f, params, upwinding)
	γ, Δx, Δt, N = params
	dz_dt = [zeros(N+1),zeros(N+1),zeros(N+1)]
	dzl_dt = [zeros(N+1),zeros(N+1),zeros(N+1)]
	if upwinding == true
		spacial_partials = [zeros(N+1),zeros(N+1),zeros(N+1),zeros(N+1),zeros(N+1)]
	else
		spacial_partials = [zeros(N+1),zeros(N+1),zeros(N+1),zeros(N+1)]
	end
	Ax = As.(x)
	derivative!(Ax, spacial_partials[4], Δx)

	boundary_coeffs = 1/2Δx*[-3.,4.,-1.]

	z_last = copy(z₀)
	z_now = copy(z₀)

	diffSys!(z_now, x, dz_dt, spacial_partials, Ax, boundary_coeffs, params, upwinding)
	z_now = z_now .+ Δt*dz_dt
	t_f = t+1001* Δt

	while t < t_f
		diffSys!(z_now, x, dz_dt, spacial_partials, Ax, boundary_coeffs, params, upwinding)
		diffSys!(z_last, x, dzl_dt, spacial_partials, Ax, boundary_coeffs, params, upwinding)
		z_last = z_now
		z_now = z_now + 3./2.*Δt*dz_dt - 1./2.*Δt*dzl_dt

		t += Δt
	end
	println(dz_dt[2])

	return z_now
end

function diffSys!(z, x, dz_dt, spacial_partials, Ax, deriv_coeffs, params, upwinding)
	ρ, v, T = z
	∂ρₜ, ∂vₜ, ∂Tₜ = dz_dt
	γ, Δx, Δt, N = params
	if upwinding == true
		∂ρₓ, ∂vₓ, ∂Tₓ, ∂Aₓ, ∂vₓ_upw = spacial_partials
		derivative_upw!(v, ∂vₓ_upw, Δx)
	else
		∂ρₓ, ∂vₓ, ∂Tₓ, ∂Aₓ = spacial_partials
	end

	derivative!(ρ, ∂ρₓ, Δx)
	derivative!(v, ∂vₓ, Δx)
	derivative!(T, ∂Tₓ, Δx)

	∂ρₜ[2:end-1] .= @view (-1. ./Ax .* (∂ρₓ.*Ax.*v .+ ρ.*∂Aₓ.*v .+ ρ.*Ax.*∂vₓ))[2:end-1]
	if upwinding == true
		∂vₜ[2:end-1] .= @view (-v .* ∂vₓ_upw .- 1./(γ*ρ) .* (∂Tₓ.*ρ .+ T.*∂ρₓ))[2:end-1] # upwinding change
	else
		∂vₜ[2:end-1] .= @view (-v .* ∂vₓ .- 1./(γ*ρ) .* (∂Tₓ.*ρ .+ T.*∂ρₓ))[2:end-1] # upwinding change
	end
	∂Tₜ[2:end-1] .= @view (-v.*∂Tₓ .- (γ-1)*(T .* ∂vₓ .+ T.*v.*∂Aₓ./Ax))[2:end-1]
	# boundary points
	∂vₜ[1] = -v[1] * sum(deriv_coeffs.*(@view v[1:3])) - 1./(γ*ρ[1]) * sum(deriv_coeffs.*((@view T[1:3]).*(@view ρ[1:3])))
	∂vₜ[end] = -v[end] * sum(-deriv_coeffs.*(@view v[end:-1:end-2])) - 1./(γ*ρ[end]) * sum(-deriv_coeffs.*((@view T[end:-1:end-2]).*(@view ρ[end:-1:end-2])))
	∂Tₜ[end] = -v[end] * sum(-deriv_coeffs.*(@view T[end:-1:end-2])) - (γ-1)*(T[end] * sum(-deriv_coeffs.*(@view v[end:-1:end-2])) + T[end]*v[end]*sum(-deriv_coeffs.*log.(As.(@view x[end:-1:end-2]))))
end

function derivative_upw!(f,∂fₓ,h)
	for i=2:length(f)-1
		if f[i] >= 0.
			∂fₓ[i] = (f[i]-f[i-1])/h
		else
			∂fₓ[i] = (f[i+1]-f[i])/h
		end
	end
end

function derivative!(f,∂fₓ,h)
	for i=2:length(f)-1
		∂fₓ[i] = (f[i+1]-f[i-1])/2h
	end
end

N = 120
Δt = 1e-4
ρ_e = 0.999
γ = 1.4
t_f = 1.
t = 0.
x = linspace(0.,1.,N+1)
Δx = 1./N

ρ = fill(ρ_e, N+1)
ρ[1] = 1.
v = fill(0., N+1)
T = fill(1., N+1)
z₀ = [ρ,v,T]
params = [γ, Δx, Δt, N]

z = AB2(z₀, x, t, t_f, params, false)
#println(z[1])
#println(z[2])
#println(z[3])

fig,(ax1,ax2,ax3) = plt[:subplots](3,1)

ax1[:plot](x, z[1], label=L"$\rho$", c="#0256dd")
ax2[:plot](x, z[2], label=L"$v$", c="#3c8409")
ax3[:plot](x, z[3], label=L"$T$", c="#f2071a")

ax1[:legend](); ax2[:legend](); ax3[:legend]()

plt[:show]()
