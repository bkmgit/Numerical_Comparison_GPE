function Reinit_M( Simplexes::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1} , M::SparseMatrixCSC{Float64,Int64}, N_of_Simplexes::Int64)

    for i = 1: N_of_Simplexes

        Simplex = Simplexes[:,i];
        idx = Mesh2Space[Simplex];

 
        

            for i = 1:3
                for j = 1:3
                    if( (idx[i]!=0) & (idx[j]!=0) )  #Om ej kant
                         M[idx[i],idx[j]] = eps(1/10000);       

                     end
                end

            end

     #   end

   

    end

        return M   

end

