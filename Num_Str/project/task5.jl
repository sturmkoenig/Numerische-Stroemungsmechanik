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
	Δt = 1e-3
	Δx = 1/(N)
	x = linspace(0,1,N+1)
	global A = Area.(x)
	global ∂A = ∂Area.(x)
	
	
	function f_pred(x)
		g = ∂t(x)
		return x + Δt*g, g 
	end
	
	function ∂f_corr(x)
		return ∂t(x, true)
	end
	
	function MC(x)
		# y is the predicted parameter
		y, ∂y= f_pred(x)
		g = x + Δt/2 * (∂y + ∂f_corr(y))
	end
	
	
	function ∂x(f, corr=false)
		g = similar(f)
	
		for i in 2:(size(f,1)-1) 
			if corr
				g[i] = (f[i+1]-f[i])/Δx
			else
				g[i] = (f[i]-f[i-1])/Δx
			end
		end
	
		g[1] = -(f[3]-4f[2]+3f[1])/2Δx
		g[end] = (f[end-2] - 4f[end-1] + 3f[end])/2Δx
		return g
	end
	
	function ∂t(x::variables, corr=false)
		ρ = ∂ρ(x.ρ, x.v, corr)
		v = ∂v(x.ρ, x.v, x.T, corr)
		T = ∂T(x.v, x.T, corr)
		erg = variables(ρ,v,T)
	end
	
	function ∂ρ(ρ,v, corr=false)
		global A
		-(1./A) .* ∂x(ρ  .* v .* A, corr)
	end
	
	function ∂v(ρ,v,T, corr=false)
		-1./(γ.*ρ) .* ∂x(ρ.*T, corr) .- v .* ∂x(v, corr)
	end
	
	function ∂T(v,T, corr=false)
		global A
		global ∂A
		-1 .* v .* ∂x(T, corr) .- (γ-1) .* (T.*∂x(v, corr) .+ T .* v .* ∂A./A)
	end

	
	function BC!(x::variables)
		x.ρ[1]=1.
		x.ρ[end]=ρe
		x.T[1]=1.
		return
	end
	
	var = variables(fill(ρe,N+1), zeros(N+1), ones(N+1))
	var.ρ[1]=1
	
	steps = 0:Δt:400
	solution = Array{variables}(size(steps,1))
	
	solution[1] = var
	
	plot(x, solution[1].v)
	for i in 2:(size(steps,1))
		solution[i] = MC(solution[i-1]) 
		BC!(solution[i])
		if (i%100) == 0
			cla()
			ylim(-0.001,0.1)
			title(i)
			plot(x,solution[i].v)	
			sleep(0.01)
		end
	end
end

integrator()
