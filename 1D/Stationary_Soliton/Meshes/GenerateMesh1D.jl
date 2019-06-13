
for N =51200
	X_L = -20; X_R = 20;

	I = X_R-X_L;
	#Nodes = [X_L+i*(I/N) for i = 0:N];
	Nodes = Array(LinRange(X_L,X_R,N+1));

	Simplices = zeros(2,N);
	for i = 1:N
	    Simplices[1,i] = i; Simplices[2,i] = i+1;
	end


	save("1dMesh_"*string(N)*".jld","Nodes",Nodes,"Simplices",Simplices)
end
