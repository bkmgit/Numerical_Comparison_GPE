function AssembleMatrix_A( Simplexes::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1},A::SparseMatrixCSC{Float64,Int64}, N_of_Simplexes::Int64 )
#Calculates du,dv
   
    for i = 1: N_of_Simplexes

        simplex = Simplexes[:,i];


        Vertexes = Nodes[:,simplex];
        J = hcat(Vertexes[:,2]-Vertexes[:,1], Vertexes[:,3]-Vertexes[:,1]);


        DetJ = det(J);

        Jinv = 1/DetJ*[J[2,2] -J[2,1]; -J[1,2] J[1,1]];

        #---Calculate local matrix-------
        Vec= hcat(Jinv*[-1.0;-1.0],Jinv*[1.0;0.0],Jinv*[0.0;1.0]);
        loc = Vec'*Vec*abs(DetJ)/2.0;


         #---Local matrix calculated------


         #Add to Global Matrix

        idx = Mesh2Space[simplex];
            for i = 1:3
                for j = 1:3       
           
                     if( (idx[i] != 0) & (idx[j] != 0) )  #Om ej kant
                        A[idx[i],idx[j]] =  A[idx[i],idx[j]]+loc[i,j];    
                     end
                end

            end

         #Added   

    end


        return A;

end




