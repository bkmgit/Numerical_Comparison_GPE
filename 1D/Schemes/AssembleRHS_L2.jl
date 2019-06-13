function AssembleRHS_L2(Simplexes::Array{Int64,2},Nodes::Array{Float64,1},Mesh2Space::Array{Int64,1} , N_of_Simplexes::Int64,Quad::Dict{String,Array{Float64,1}},b,fuck)
#Calculates <|phi|^2 phi_j, phi_i>
    
    N_gp = Int(Quad["N_gp"][1])
#Help functions

 f(X::Float64,C1::Complex{Float64},C2::Complex{Float64}) = (1-X)*C1 +C2*X         
#---------------------------

#initilise variables -------
    phi = [0.0 0.0];
    U = 0.0+0.0im;


    
    loc = zeros(2,1)*0im;


    for i = 1: N_of_Simplexes

        Simplex = Simplexes[:,i];

        Vertexes = Nodes[Simplex];
         
        J = Vertexes[2]-Vertexes[1];


        DetJ = abs(J);


        

        fill!(loc,0);

        for gp = 1:N_gp
            X = Quad["X"][gp]; ω =Quad["ω"][gp];

            phi = Φ(X,phi);
            x = J*X+Vertexes[1];
  	   for i = 1:2
                 loc[i] +=  fuck(x)*phi[i]*ω;
            end
        end
            loc = loc*DetJ;

    
         #Add to Global Matrix
         #Add to Global Matrix
         idx = Mesh2Space[Simplex];


            for i = 1:2
                if( (idx[i]!=0)  )  #Om ej kant
                         b[idx[i]] = b[idx[i]]+loc[i];       

                end

            end


    end



        return b

end

