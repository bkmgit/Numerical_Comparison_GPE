function Define_Space(BoundaryX::Float64, BoundaryY::Float64,N::Int64,Nodes::Array{Float64,2})
    Mesh2Space=zeros(Int64,N);


   


    index = 1;
    for i = 1:N
        node = Nodes[:,i];
        if( (node[1] == -BoundaryX) | (node[2] ==-BoundaryY) )

        elseif( (node[1] == BoundaryX) | (node[2] == BoundaryY) )

        else
            Mesh2Space[i] = index;
            index = index+1;
        end
    end

    #Mesh2Space=convert(Array{Int64,2},Mesh2Space);

    SpaceSize = index-1;

    
    return Mesh2Space,SpaceSize
    
end
