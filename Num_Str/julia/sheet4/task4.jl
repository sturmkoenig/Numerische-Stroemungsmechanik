function dT/dt(T, x, Q, h)
	return central_difference_2od(T, x, h) + Q(x)
end

function euler_explicit(x₀, f)
    x₁ = x₀ + Δt*f(x₀)
end

function Q(x)
	sinpi(x)
end

function forward_difference(f,x,h)
	return f(x+h)-f(x)	
end

function central_difference_2od(f, x, h)
	return (f(x+h) -2f(x)+f(x-h))/h^2
end

function FTCS(T, x, Q, h)
	return T+dt* dT/dt(T, x, Q, h)
end


for i in 1:length(steps)
	save[i] = central_difference_2od(f,i*h, h)
end

plot(steps, save)
	
