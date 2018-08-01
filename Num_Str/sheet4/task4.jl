function dT(T, x, h)
	return central_difference_2od(T, h) + Q(x)
end

function Q(x)
	sinpi.(x)
end

function euler(T, x, h)
	return T + h*dT(T, x, h)
end

function central_difference_2od(f, h)

	return (circshift(f, 1) + circshift(f, -1) - 2f)/h^2
end

N = 10
x = linspace(0, 1, N)
T = zeros(N)	

for i in 1:10
	T=euler(T, x, 0.01)
	println(T)
end

plot(x, T)

plt[:show]()