function AssembleMatrix_M_RE( Simplexes::Array{Int64,2},Nodes::Array{Float64,1},Mesh2Space::Array{Int64,1} ,p, M::SparseMatrixCSC{Float64,Int64}, N_of_Simplexes::Int64,Quad::Dict{String,Array{Float64,1}})
#Calculates <|phi|^2 phi_j, phi_i>
    
    N_gp = Int(Quad["N_gp"][1])
#Help functions

#initilise variables -------
    phi = [0.0 0.0];

    
    loc = zeros(2,2);

#--------------------------------
it = 1;
    for i = 1: N_of_Simplexes

        Simplex = Simplexes[:,i];

        Vertexes = Nodes[Simplex];
         
        J = Vertexes[2]-Vertexes[1];

        idx = Mesh2Space[Simplex];

        DetJ = abs(J);

        fill!(loc,0);

        for gp = 1:N_gp
            X = Quad["X"][gp]; ω =Quad["ω"][gp];
            Quad_Poly = 2*(X-0.5)*(X-1)*p[1,it]+2*X*(X-0.5)*p[2,it]+4*(1-X)*X*p[3,it];
            phi = Φ(X,phi);
 
            for i = 1:2
                for j = 1:2
                    loc[i,j] +=  Quad_Poly*phi[i]*phi[j]*ω;
                end
            end
        end
            loc = loc*DetJ;

    
         #Add to Global Matrix


            for i = 1:2
                for j = 1:2
                    if( (idx[i]!=0) & (idx[j]!=0) )  #Om ej kant
                         M[idx[i],idx[j]] =  M[idx[i],idx[j]]+loc[i,j];       

                     end
                end

            end

it+=1;
    end



        return M   

end

