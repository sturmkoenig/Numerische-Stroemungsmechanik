import Base: +
import Base: -
import Base: *

struct variables
	ρ::Vector{Float64}
	v::Vector{Float64}
	T::Vector{Float64}
end

function +(A::variables, B::variables)
	ρ = A.ρ .+ B.ρ
	v = A.v .+ B.v
	T = A.T .+ B.T
	erg = variables(ρ,v,T)
end

function -(A::variables, B::variables)
	ρ = A.ρ .- B.ρ
	v = A.v .- B.v
	T = A.T .- B.T
	erg = variables(ρ,v,T)
end

function *(A, B::variables)
	ρ = A * B.ρ
	v = A * B.v
	T = A * B.T
	erg = variables(ρ,v,T)
end




function integrator()

	function Area(x) 
		if 0 <= x <= .5
			return 1. + 2.2(x - 0.5)^2
		elseif .5 < x <= 1
			return 1. + 0.22(x - 0.5)^2
		end
	end
	
	function ∂Area(x)
		if 0<= x <= .5
			return 4.4(x-.5)
		elseif .5 < x <= 1
			return .44(x-.5)
		end
	end

	ρe = 0.999
	γ = 1.4
	N = 120
	Δt = 1e-4
	Δx = 1/(N)
	upwinding = false
	x = linspace(0,1,N+1)
	global A = Area.(x)
	global ∂A = ∂Area.(x)
	global ∂lnA = ∂A./A

	function ∂x(f, upw=false)
		g = similar(f)
	
		for i in 2:(size(f,1)-1) 
			if upw
				if f[i]<0
					g[i] = (f[i+1]-f[i])/Δx
				else
					g[i] = (f[i]-f[i-1])/Δx
				end
			else
				g[i] = (f[i+1]-f[i-1])/(2Δx)
			end
		end
	
		g[1] = -(f[3]-4f[2]+3f[1])/2Δx
		g[end] = (f[end-2] - 4f[end-1] + 3f[end])/(2Δx)
		return g
	end
	
	function ∂t(x::variables)
		ρ = -(1./A) .* ∂x(x.ρ  .* x.v .* A)
		v = -1./(γ.*x.ρ) .* ∂x(x.ρ.*x.ρ) .- x.v .* ∂x(x.v, true)
		T = -1 .* x.v .* ∂x(x.T) .- (γ-1) .* (x.T.*∂x(x.v) .+ x.T .* x.v .* ∂lnA)
		erg = variables(ρ,v,T)
	end
		
	function euler_explicit(x₀)
	    x₁ = x₀ .+ Δt.*∂t(x₀)
	end
	
	function AB1(x₀, x₁)
	    x₂ = x₁ .+ 1/2*Δt.*(3∂t(x₁) .- ∂t(x₀))
	end
	
	
	function BC!(x::variables)
		x.ρ[1]=1.
		x.ρ[end]=ρe
		x.T[1]=1.
		return
	end
	
	var = variables(fill(ρe,N+1), zeros(N+1), ones(N+1))
	var.ρ[1]=1
	
	steps = 0:Δt:1
	solution = Array{variables}(size(steps,1))
	
	solution[1] = var
	solution[2] = euler_explicit(solution[1])
	
	fig = figure()
	ax = axes()
	
	BC!(solution[1])
	plot(x, solution[1].v)
	for i in 3:(size(steps,1))
		solution[i] = AB1(solution[i-2], solution[i-1]) 
		BC!(solution[i])
		if (i%100) == 0
			cla()
			ylim(-0.001,0.01)
			title(i)
			plot(x,solution[i].v)	
			sleep(0.01)
		end
	end
end

integrator()
