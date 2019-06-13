function L2_Projection(Simplexes::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1} ,N_of_Simplexes::Int64,N_of_Nodes::Int64,SpaceSize::Int64,Quad::Dict{String,Array{Float64,1}},sitp,Boundary)
#Calculates <|phi|^2 phi_j, phi_i>
    
    N_gp = Int(Quad["N_gp"][1])
#Help functions

 
#initilise variables -------
    phi = [0.0 0.0 0.0]; 
    U = 0.0+0.0im;


    loc = zeros(Complex{Float64},3,3);#zeros(3,1)+0im;

#--------------------------------


p = zeros(Complex{Float64},SpaceSize,1);#zeros(SpaceSize,1)+0im;

    for i = 1: N_of_Simplexes

        Simplex = Simplexes[:,i];

        Vertexes = Nodes[:,Simplex];
         
        J = hcat(Vertexes[:,2]-Vertexes[:,1], Vertexes[:,3]-Vertexes[:,1]);

        DetJ = det(J);
        
        idx = Mesh2Space[Simplex];



        fill!(loc,0);

        for gp = 1:N_gp
            X = Quad["X"][gp]; ω =Quad["ω"][gp];Y = Quad["Y"][gp];

            phi = Φ(X,Y,phi);
            x = Vertexes[:,1]+J*[X;Y];
	    if( (abs(x[1])<Boundary) && (abs(x[2])<Boundary) ) w = sitp(x[1],x[2]); else w = 0.0; end
            for i = 1:3
                loc[i]+= w*ω*phi[i];
            end
        end
            loc = loc*abs(DetJ)/2.0;

    
         #Add to Global Matrix
         #Add to Global Matrix
 

            for i = 1:3
                if( (idx[i]!=0)  )  #Om ej kant
                         p[idx[i]] = p[idx[i]]+loc[i];       

                end

            end

    end



        return p

end

