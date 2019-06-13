function AssembleMatrix_M_RE( Simplexes::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1} ,ϕ::Array{Float64,2}, M::SparseMatrixCSC{Float64,Int64}, N_of_Simplexes::Int64,Quad::Dict{String,Array{Float64,1}})
#Calculates <|phi|^2 phi_j, phi_i>
    N_gp = Int(Quad["N_gp"][1])
#Help functions

#initilise variables -------
    phi = [0.0 0.0 0.0];
    u = 0.0+0.0im;

    N = size(Nodes,2)+0.0im;

    loc = zeros(3,3);


    for i = 1: N_of_Simplexes

        Simplex = Simplexes[:,i];
        idx = Mesh2Space[Simplex];

        Vertexes = Nodes[:,Simplex];
         
        J = hcat(Vertexes[:,2]-Vertexes[:,1], Vertexes[:,3]-Vertexes[:,1]);

        #DetJ = abs(det(J));
 	DetJ = det(J);
        #---Calculate local matrix-------
                 

        

        
        fill!(loc,0);
        
        C = ϕ[:,i];
  # if(sum(abs.(C))>10.0^-16)     
        for gp = 1:N_gp
            X = Quad["X"][gp]; Y = Quad["Y"][gp];ω =Quad["ω"][gp];

            phi = Φ(X,Y,phi);
          
            nodal_basis = [phi[1]*(2*phi[1]-1),4*phi[1]*phi[2],phi[2]*(2*phi[2]-1),4*phi[1]*phi[3] ,4*phi[2]*phi[3],phi[3]*(2*phi[3]-1)]';
            
            for i = 1:3
                for j = 1:3
                    loc[i,j] +=(nodal_basis*C)*phi[i]*phi[j]*ω 
                end
            end
            
         end
        loc = loc*abs(DetJ)/2.0


        

            for i = 1:3
                for j = 1:3
                    if( (idx[i]!=0) & (idx[j]!=0) )  #Om ej kant
                         M[idx[i],idx[j]] +=  loc[i,j]#+eps(1/10000);       

                     end
                end

            end

     #   end

   

    end

        return M   

end

