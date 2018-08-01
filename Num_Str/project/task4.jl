import Base: +
import Base: -
import Base: *

ρe = 0.999
γ = 1.4
N = 320
Δt = 1e-4
Δx = 1/(N)

struct variables
	ρ
	v
	T
end

function Area(x) 
	if 0 <= x <= .5
		return 1 + 2.2(x - 0.5)^2
	elseif .5 < x <= 1
		return 1 + 0.22(x - 0.5)^2
	end
end

x = linspace(0,1,N+1)
A = Area.(x)

function ∂x(f)
	g = similar(f)

	for i in 2:(size(f,1)-1)
		g[i] = (f[i+1]-f[i-1])/2Δx
	end

	g[1] = -(f[3]-4f[2]+3f[1])/2Δx
	g[end] = -(f[end-2] - 4f[end-1] + 3f[end])/2Δx
	return g
end

function ∂ρ(ρ,v,A)
	-(1./A) .* ∂x(ρ  .* v .* A)
end

function ∂v(ρ,v,T)
	-1./(γ.*ρ) .* ∂x(ρ.*T) .- v .* ∂x(v)
end

function ∂T(v,T,A)
	-1 .* v .* ∂x(T) .- (γ-1) .* (T.*∂x(v) .+ T .* v .* A./∂x(A))
end


function euler_explicit(x₀)
    x₁ = x₀ .+ Δt.*∂t(x₀)
end

function AB1(x₀, x₁)
    x₂ = x₁ .+ 1/2*Δt.*(3∂t(x₁) .- ∂t(x₀))
end

function AB2(x₀, x₁, x₂)
    x₃ = x₂ .+ 1/12*Δt*(23∂t(x₂) .- 16∂t(x₁) .+ 5∂t(x₀))
end

function ∂t(x::variables)
	ρ = ∂ρ(x.ρ, x.v, A)
	v = ∂v(x.ρ, x.v, x.T)
	T = ∂T(x.v, x.T, A)
	erg = variables(ρ,v,T)
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

function BC!(x::variables)
	x.ρ[1]=1.
	x.ρ[end]=ρe
	x.T[1]=1.
	return
end

var = variables(fill(ρe,N+1), zeros(N+1), ones(N+1))
var.ρ[1]=1

steps = 0:Δt:1
global solution = variables[]

push!(solution, var)
push!(solution, euler_explicit(solution[1]))

fig = figure()
ax = axes()
plot(x,A)

BC!(solution[end])

#plot(x, solution[end].v)
#for i in 2:5400
#	if (i%100) == 0
#		cla()
#		ylim(-0.01,0.005)
#		title(i)
#		plot(x, solution[end].v)
#		sleep(0.01)
#	end
#	push!(solution, AB1(solution[i-1], solution[i])) 
#	BC!(solution[end])
#end

println(solution[end].ρ)
println(solution[end].v)
println(solution[end].T)

#plot(x, solution[end].ρ)
