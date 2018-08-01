function dT(T, x, N, t)
	Ret = central_difference_2od(T, 1/N)
	Ret[1]=0
	Ret[N]=sinpi.(4/3*50t)
	return Ret
end

function euler(T, x, dt, N, t)
	return T + dt*dT(T, x, N, t)
end

function central_difference_2od(f, h)
	return (circshift(f, 1) - 2f + circshift(f, -1) )/h^2
end

N = 100
x = linspace(0, 1, N)
T = zeros(N)
dt = 2e-5
time = collect(0:dt:2)

for t in time
	T=euler(T, x, dt,  N, t)
end

plot(x, T)

plt[:show]
