function dT(T, x, N, t)
	Ret = D * central_difference_2od(T, 1/N) + v*central_difference_1od(T, 1/N)
end

function FCTS(T, x, dt, N, t)
	T = T + dt*dT(T, x, N, t)
	T[1]=0
	T[N]=1
	return T
end

function central_difference_1od(f, h)
	return (circshift(f, 1) - circshift(f, -1))/h
end

function central_difference_2od(f, h)
	return (circshift(f, 1) - 2f + circshift(f, -1) )/h^2
end

N = 1000
dt = 0.0001
v=1
D=1
x = linspace(0, 1, N)
T = linspace(0, 1, N)

steps = collect(0:dt:1)

for t in steps
	T=FCTS(T, x, dt, N, t)
end

print(T)
plot(x, T)
plt[:show]()
