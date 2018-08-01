function d(x)
	return x^2
end

function g(x)
	return x^3
end

function forward_difference(f,x,h)
	return f(x+h)-f(x)	
end

function central_difference_2od(f, x, h)
	return (f(x+h) -2f(x)+f(x-h))/h^2
end


h = 0.3
p = [0]
space = 10
steps = linspace(0, 10, 10/0.1)
save = similar(steps)

for i in 1:length(steps)
	save[i] = central_difference_2od(f,i*h, h)
end

plot(steps, save)
	
