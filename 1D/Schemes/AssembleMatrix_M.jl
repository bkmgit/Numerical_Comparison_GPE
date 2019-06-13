function  AssembleMatrix_M(  Simplexes::Array{Int64,2},Nodes::Array{Float64,1},Mesh2Space::Array{Int64,1} ,weight::Function,M::SparseMatrixCSC{Float64,Int64},N_of_Simplexes::Int64,Quad::Dict{String,Array{Float64,1}})
#Calculates <weight[x]* phi_j, phi_i>
  
    N_gp = Int(Quad["N_gp"][1])
#Help functions

    
 #initialise stuff   
    loc = zeros(2,2);
    phi= [0.0 0.0];

#---------------------

    for i = 1:N_of_Simplexes
    
        Simplex = Simplexes[:,i];
   

        Vertexes = Nodes[Simplex];
        J = Vertexes[2]-Vertexes[1];


        DetJ = abs(J);
        F(X) = J*X+Vertexes[1];
 

        
               

        fill!(loc,0);
        for gp = 1:N_gp
            X = Quad["X"][gp]; ω =Quad["ω"][gp];

            phi = Φ(X,phi);
            w = weight(F( X) )
            for i = 1:2
                for j = 1:2
                    loc[i,j] = loc[i,j]+ w*phi[i]*phi[j]*ω;
                end
            end
        end
        loc = loc*DetJ
        
     #---Local matrix calculated------
       
     #Add to Global Matrix
  
        idx = Mesh2Space[Simplex];
        
        for i = 1:2
                
            for j = 1:2
                
                 if( (idx[i]!=0) & (idx[j]!=0) ) #Om ej kant
                     M[idx[i],idx[j]] =  M[idx[i],idx[j]]+loc[i,j];                
                 end
            end
            
        end
      
     #Added   

 
end

    return M;
end


