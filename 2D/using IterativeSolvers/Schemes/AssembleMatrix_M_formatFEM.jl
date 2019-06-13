function AssembleMatrix_M_formatFEM( Simplexes::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1} ,weight::Array{Complex{Float64},1}, M::SparseMatrixCSC{Float64,Int64}, N_of_Simplexes::Int64,Quad::Dict{String,Array{Float64,1}})
#Calculates <|phi|^2 phi_j, phi_i>
    
    N_gp = Int(Quad["N_gp"][1])
#Help functions

   
#---------------------------

#initilise variables -------
    phi = [0.0 0.0 0.0];

    N = size(Nodes,2);

    C = zeros(Complex{Float64},3)

    loc = zeros(Complex{Float64},3,3)

#--------------------------------


    for i = 1: N_of_Simplexes

        simplex = Simplexes[:,i];

        Corners = Nodes[:,simplex];
        J = hcat(Corners[:,2]-Corners[:,1], Corners[:,3]-Corners[:,1]);

        #DetJ = abs(det(J));
        DetJ = det(J);
        idx = Mesh2Space[simplex];
        
        max = 0.0;

        fill!(C,0);
        for i = 1:3 if(idx[i]!=0) C[i] = weight[idx[i]]; if(max<abs(C[i])) max = abs(C[i]); end; end; end;


  #      if(max>10.0^-16);
            fill!(loc,0);

            for gp = 1:N_gp
            X = Quad["X"][gp]; Y = Quad["Y"][gp];ω =Quad["ω"][gp];

                phi = Φ(X,Y,phi);
                w = abs(U_h(X,Y,C))^2;
                for i = 1:3
                    for j = 1:3
                        loc[i,j] +=  w*phi[i]*phi[j]*ω;
                    end
                end
            end
                loc = loc*abs(DetJ)/2.0;


             #Add to Global Matrix
             #Add to Global Matrix


                for i = 1:3
                    for j = 1:3
                        if( (idx[i]!=0) & (idx[j]!=0) )  #Om ej kant
                             M[idx[i],idx[j]] =  M[idx[i],idx[j]]+loc[i,j];       

                         end
                    end

                end
     #       end

    end



        return M   

end

