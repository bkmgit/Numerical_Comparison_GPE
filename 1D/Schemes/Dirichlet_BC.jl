function Dirichlet_BC(X_L::Float64, X_R::Float64,N::Int64,Nodes::Array{Float64,1})
    
    
    Mesh2Space=zeros(Int64,N);
    fill!(Mesh2Space,0);


    index = 1;
    
    # NO DIRICHLET BC!
  
    for i = 1:N
        node = Nodes[i];
        if( (node == X_L) | (node ==X_R) )

        else
            Mesh2Space[i] = index;
            index = index+1;
        end
    end

    SpaceSize = index-1;
 #Mesh2Space=convert(Array{Int64,2},Mesh2Space);
 #=
    for i = 1:N
        Mesh2Space[i] = index;
        index= index+1;
    end
    SpaceSize = index-1;
=#
 
    return Mesh2Space,SpaceSize
    
end
