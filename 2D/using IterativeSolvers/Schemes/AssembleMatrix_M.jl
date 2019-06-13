function  AssembleMatrix_M(  Simplexes::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1} ,weight::Function,M::SparseMatrixCSC{Float64,Int64},N_of_Simplexes::Int64,Quad::Dict{String,Array{Float64,1}})
#Calculates <weight[x]* phi_j, phi_i>

N_gp = Int(Quad["N_gp"][1])


 #initialise stuff   
    loc = zeros(3,3);
    phi= [0.0 0.0 0.0];

    for i = 1: N_of_Simplexes
	

        simplex = Simplexes[:,i];
	
        Vertexes = Nodes[:,simplex];
        J = hcat(Vertexes[:,2]-Vertexes[:,1], Vertexes[:,3]-Vertexes[:,1]);


        DetJ = det(J);

        #---Calculate local matrix-------


        loc*=0;		       

		for gp = 1:N_gp
		    X = Quad["X"][gp]; Y = Quad["Y"][gp];ω =Quad["ω"][gp];

		    phi = Φ(X,Y,phi);

		    w = weight(J*[X;Y]+Vertexes[:,1])
		    for i = 1:3
		        for j = 1:3
		            loc[i,j] += w*phi[i]*phi[j]*ω;
		        end
		    end
		end
		loc *= abs(DetJ)/2.0
  
     #---Local matrix calculated------

	idx = Mesh2Space[simplex];
     
       
	#Add to Global Matrix

        
            for i = 1:3
                for j = 1:3       
           
                     if( (idx[i] != 0) & (idx[j] != 0) )  #Om ej kant 
                        M[idx[i],idx[j]] += loc[i,j];    
                     end
                end

            end

     #Added   

 
end

    return M;
end

