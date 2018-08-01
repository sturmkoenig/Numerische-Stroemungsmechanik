function f(x)
	return sin.(5x) 
end

function forward_difference(f,x,h)
	return f(x+h)-f(x)	
end

function central_difference_2od(f, x, h)
	return (f(x+h) -2f(x)+f(x-h))/h^2
end

I = linspace(0, 4pi, 100)
save = similar(I)

for i in 1:length(I)
	save[i] = central_difference_2od(f,I[i], 4pi/100)
end

plot(I, save)
plot(I, -25*f(I))
	
