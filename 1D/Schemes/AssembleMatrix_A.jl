function AssembleMatrix_A( Simplexes::Array{Int64,2},Nodes::Array{Float64,1},Mesh2Space::Array{Int64,1},A::SparseMatrixCSC{Float64,Int64}, N_of_Simplexes::Int64 )
#Calculates du,dv
    
    
    

    for i = 1: N_of_Simplexes

        Simplex = Simplexes[:,i];


        Vertexes = Nodes[Simplex];
        
        J = Vertexes[2]-Vertexes[1];


        DetJ = abs(J);

        Jinv = 1/J;

        #---Calculate local matrix-------
        Vec = Jinv*[-1.0 1.0];
        loc = Vec'*Vec*DetJ;


         #---Local matrix calculated------


         #Add to Global Matrix
        idx = Mesh2Space[Simplex];
            for i = 1:2
                for j = 1:2
                     if( (idx[i]!=0) & (idx[j]!=0) )  #Om ej kant
                         A[idx[i],idx[j]] =  A[idx[i],idx[j]]+loc[i,j];                
                     end
                end

            end

         #Added   

    end


        return A;

end




