function central_difference_2od(f, h)
	return (circshift(f, 1) - 2f + circshift(f, -1) )/h^2
end

function dT(T, r, N, t)
	Ret = central_difference_2od(T, 1/N)
	Ret[1]=0
	Ret[N]=sinpi.(4/3*50t)
	return Ret
end

Tw=0
T0=1
r0=1
