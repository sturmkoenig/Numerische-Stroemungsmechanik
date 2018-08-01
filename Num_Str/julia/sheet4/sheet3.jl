using PyPlot
using FastAnonymous

Δt = 0.001
ε = 0.01
k = 1
M = 1

function adams_bashforth_1(x₀, x₁)
    x₂ = x₁ + 1/2*Δt*(3f(x₁) - f(x₀))
end

function adams_bashforth_2(x₀, x₁, x₂)
    x₃ = x₂ + 1/12*Δt*(23f(x₂) - 16f(x₁) + 5f(x₀))
end

function adams_bashforth_3(x₀, x₁, x₂, x₃)
    x₄ = x₃ + 1/24*Δt*(55f(x₃) - 59f(x₂) + 37f(x₁) - 9f(x₀))
end

function crank_nicolson(x₀, ε, f, df)
    g = @anon x -> x - x₀ - 1/2*Δt*(f(x) + f(x₀))
    dg = @anon x -> eye(length(x₀), length(x₀)) - 1/2*Δt*df(x)
    x₁ = newtons_method(x₀, ε, g, dg)
end

function euler_implicit(x₀, ε, f, df)
    g = @anon x-> x - Δt*f(x) - x₀
    dg = @anon x-> eye(length(x₀), length(x₀)) - Δt*df(x)
    return newtons_method(x₀, ε, g, dg) 
end


function euler_explicit(x₀, f)
    x₁ = x₀ + Δt*f(x₀)
end

function velocity_verlet(x₀, v₀, f)
    x₁ = x₀ + Δt*(v₀ + 1/2*Δt*f(x₀))
    v₁ = (x₁-x₀)/Δt + 1/2*Δt*f(x₁)
    return x₁, v₁ 
end

function newtons_method(x₀, ε, f, df)
    x₁ = x₀ + inv(df(x₀))*f(x₀)
    while(norm(x₁-x₀) > ε)
        x₀ = x₁
        x₁ = x₀ + inv(df(x₀))*f(x₀)
    end
    return x₁
end


function f(x)
    return [x[2], -k/M*x[1]]
end

function df(x)
    return [0 -k/M; 1 0]
end

function g(x)
    return -k/M*x 
end

fig = figure()
t = 100000
x = zeros(2, t)
x[:, 1] = [0,1]
x[:, 2] = euler_explicit(x[:,1], f)
x[:, 3] = adams_bashforth_1(x[:,1], x[:,2])
for i in 4:t 
    x[:, i] = adams_bashforth_2(x[:,i-3], x[:,i-2], x[:,i-1])
    println(x[:,i])
end

print(x[1,:])
#plot(1:t, x[1,:])

plt[:show]()
