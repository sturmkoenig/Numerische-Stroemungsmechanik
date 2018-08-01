function dT(T,Nx, Ny)
	 return central_difference_2od(T, 1/Nx, 1) + central_difference_2od(T, 1/Ny, 2) - Pe * (vx * central_difference_1od(T, 1/Nx, 1) + vy * central_difference_1od(T, 1/Ny, 2) ) 
end

function FCTS(T, Nx, Ny)
	T = T + dt*dT(T, Nx, Ny)
	boundary_condition(T)
	return T
end

function boundary_condition(T)
	for i in 1:size(T, 1)
		T[i, 1] = T[i, size(T,2)] = 0
	end
	for i in 1:size(T, 2) 
		T[1, i] = T[size(T, 1),i] = 1 	
	end
end

function central_difference_1od(f, h, dim)
	return (circshift(f, dim) - circshift(f, -dim))/h
end

function central_difference_2od(f, h, dim)
	return (circshift(f, dim) - 2f + circshift(f, -dim) )/h^2
end

L = 1
dt = 0.001
Nx = 50
Ny = 50
Pe = 2
vx = 1
vy = 2

Y = linspace(0, L, Ny)

T = zeros(Nx, Ny)
time = collect(0:dt:2)

for t in time
	FCTS(T, Nx, Ny)
end

imshow(T)
plt[:show]()
