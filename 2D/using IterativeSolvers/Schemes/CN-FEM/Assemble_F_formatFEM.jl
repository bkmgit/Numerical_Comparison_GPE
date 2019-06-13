function Assemble_F_formatFEM( Simplices::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1} ,U_it::Array{Complex{Float64},1},U0::Array{Complex{Float64},1}, F_real::Array{Float64,1},F_imag::Array{Float64,1}, N_of_Simplices::Int64,Quad::Dict{String,Array{Float64,1}})
#Calculates <|phi|^2 phi_j, phi_i>
    
    
    N_gp = Int(Quad["N_gp"][1])

    phi = [0.0 0.0 0.0];
    u = 0.0+0.0im;
    u0 = 0.0+0.0im;

    loc = zeros(3)*0im;
    
    C = zeros(3)*0im;
    C_0 = zeros(3)*0im;

    for i = 1: N_of_Simplices

        simplex = Simplices[:,i];
        
        idx = Mesh2Space[simplex];
        
        Vertex = Nodes[:,simplex];
        J = hcat(Vertex[:,2]-Vertex[:,1], Vertex[:,3]-Vertex[:,1]);

        DetJ = abs(det(J));

        
        fill!(C,0);
        fill!(C_0,0);
        fill!(loc,0);
        
        for j = 1:3 if(idx[j]!=0) C[j] = U0[idx[j]]; C_0[j] = U_it[idx[j]]; end end

        for gp = 1:N_gp
            X = Quad["X"][gp]; Y = Quad["Y"][gp];ω =Quad["ω"][gp];
            phi = Φ(X,Y,phi);
            u = U_h(X,Y,C);
            u0 = U_h(X,Y,C_0);
            for i = 1:3
                loc[i] += (abs(u)^2+abs(u0)^2)*(u+u0)*phi[i]*ω
            end
        end
            loc = loc*DetJ/2.0;

    
         #Add to Global Matrix
         #Add to Global Matrix
         idx = Mesh2Space[simplex];


            for i = 1:3
                    if( (idx[i]!=0) )  #Om ej kant
                         F_real[idx[i]] += real(loc[i]);
                         F_imag[idx[i]] += imag(loc[i]);  
                    end

            end


    end


        
        return F_real, F_imag;

end

