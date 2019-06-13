function Assemble_F_formatFEM_theta( Triangles::Array{Int64,2},Nodes::Array{Float64,1},Mesh2Space::Array{Int64,1} ,U_it::Array{Complex{Float64},1},U0::Array{Complex{Float64},1}, F_real::Array{Float64,1},F_imag::Array{Float64,1}, N_of_Triangles::Int64,θ::Float64,Quad::Dict{String,Array{Float64,1}})
#Calculates <|phi|^2 phi_j, phi_i>
    
    
    N_gp = Int(Quad["N_gp"][1])

    
 U(X::Float64,C1::Complex{Float64},C2::Complex{Float64}) = (1-X)*C1 +C2*X 

    phi = [0.0 0.0];
    u = 0.0+0.0im;
    u0 = 0.0+0.0im;

#THIS IS STUPID, REPLACE SOMETIME
    N = length(Nodes);
    #translate solution points to whole mesh 
    Data = zeros(Complex{Float64},N);
    
    Data_0 =zeros(Complex{Float64},N);
    
    
    loc = zeros(Complex{Float64},2);

    for i = 1:N
       if( Mesh2Space[i] != 0)
          Data[i] = U_it[Mesh2Space[i]]; 
          Data_0[i] = U0[Mesh2Space[i]];
       end
    end


    for i = 1: N_of_Triangles

        triangle = Triangles[:,i];

        Vertex = Nodes[triangle];
        J =Vertex[2]-Vertex[1];

        DetJ = abs(J);



        C1 = Data[triangle[1]];
        C2 = Data[triangle[2]];


        C1_0 = Data_0[triangle[1]];
        C2_0 = Data_0[triangle[2]];



        fill!(loc,0);


        for gp = 1:N_gp
            X = Quad["X"][gp]; ω =Quad["ω"][gp];
            phi = Φ(X,phi);
            u = U(X,C1,C2);
            u0 = U(X,C1_0,C2_0);
            uθ = θ*u+(1-θ)*u0;
            for i = 1:2
                loc[i] += abs(uθ)^2*uθ*phi[i]*ω
            end
        end
            loc = loc*DetJ;

    
         #Add to Global Matrix
         #Add to Global Matrix
         idx = Mesh2Space[triangle];


            for i = 1:2
                    if( (idx[i]!=0) )  #Om ej kant
                         F_real[idx[i]] += real(loc[i]);
                         F_imag[idx[i]] += imag(loc[i]);  
                    end
            end


    end
    
    return F_real, F_imag;

end

