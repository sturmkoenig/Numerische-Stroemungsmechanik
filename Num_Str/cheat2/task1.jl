length = Int64(floor(pi/0.1))
M = zeros(length)
N = zeros(length)
for i in 1:length
	M[i]=i*0.1
	N[i]=sin(i*0.1)
end

for i in 1:length
	print("$(M[i])	$(N[i])\n")
end


