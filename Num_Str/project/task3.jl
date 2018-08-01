function A(x) 
	if 0 <= x <= .5
		return 1 + 2.2(x - .5)^2
	elseif .5 <= x <= 1
		return 1 + .22(x-.5)^2
	end
end

x = linspace(0, 1, 1000)
plot(x, A.(x))
