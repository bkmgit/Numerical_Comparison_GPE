function AssembleMatrix_J_F_formatFEM( Simplices::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1} ,U_it::Array{Complex{Float64},1},U0::Array{Complex{Float64},1},
M_1::SparseMatrixCSC{Float64,Int64},
M_2::SparseMatrixCSC{Float64,Int64},
M_3::SparseMatrixCSC{Float64,Int64},
M_4::SparseMatrixCSC{Float64,Int64}, 
N_of_Simplices::Int64,
Quad::Dict{String,Array{Float64,1}})
#Calculates <|phi|^2 phi_j, phi_i>


N_gp=( Int(Quad["N_gp"][1]))
    
function Φ(X::Float64,Y::Float64,phi::Array{Float64}) 
       
    phi[1]=1-X-Y
    phi[2]= X;
    phi[3]= Y;
    return phi
 end
    
 U(X::Float64,Y::Float64,C::Array{Complex{Float64}}) = (1-X-Y)*C[1] +C[2]*X + C[3]*Y        

    phi = [0.0 0.0 0.0];
    u = 0.0+0.0im;
    u0 = 0.0+0.0im;

    C1 = zeros(3)*0im;
    C2 = zeros(3)*0im;

    
    loc_1 = Array{Float64,2}(3,3);
    loc_2 = Array{Float64,2}(3,3);
    loc_3 = Array{Float64,2}(3,3);
    loc_4 = Array{Float64,2}(3,3);


    for i = 1: N_of_Simplices

        simplex = Simplices[:,i];

        Vertex = Nodes[:,simplex];
        J = hcat(Vertex[:,2]-Vertex[:,1], Vertex[:,3]-Vertex[:,1]);

        DetJ = abs(det(J));

	idx = Mesh2Space[simplex];
	
	fill!(C1,0)
	fill!(C2,0)
	
	for ii = 1:3 if(idx[ii]!=0) C1[ii]=U_it[idx[ii]]; C2[ii]=U0[idx[ii]]; end
            
      if(sum(abs.(C1)+abs.(C2))>=10.0^-16)


        fill!(loc_1,0);
        fill!(loc_2,0);
        fill!(loc_3,0);
        fill!(loc_4,0);

        for gp = 1:N_gp
            X = Quad["X"][gp]; Y = Quad["Y"][gp];ω =Quad["ω"][gp];

            phi = Φ(X,Y,phi);
            u = U_h(X,Y,C1);
            u0 = U_h(X,Y,C2);
            for i = 1:3
                for j = 1:3
                    loc_1[i,j] +=  (2*real(u)*(real(u)+real(u0)) + abs(u)^2+abs(u0)^2)*phi[i]*phi[j]*ω;
                
                    loc_2[i,j] += 2*imag(u)*(real(u)+real(u0))*phi[i]*phi[j]*ω;

                    loc_3[i,j] += 2*real(u)*(imag(u)+imag(u0))*phi[i]*phi[j]*ω;

                    loc_4[i,j] += (2*imag(u)*(imag(u)+imag(u0)) + abs(u)^2+abs(u0)^2)*phi[i]*phi[j]*ω;
                end
            end
        end
            loc_1 = loc_1*DetJ/2.0;
            loc_2 = loc_2*DetJ/2.0;
            loc_3 = loc_3*DetJ/2.0;
            loc_4 = loc_4*DetJ/2.0;

    
         #Add to Global Matrix
         #Add to Global Matrix
         idx = Mesh2Space[simplex];


            for i = 1:3
                for j = 1:3
                    if( (idx[i]!=0) & (idx[j]!=0) )  #Om ej kant
                         M_1[idx[i],idx[j]] +=  loc_1[i,j];  
			 M_2[idx[i],idx[j]] +=  loc_2[i,j];       
			 M_3[idx[i],idx[j]] +=  loc_3[i,j];
			 M_4[idx[i],idx[j]] +=  loc_4[i,j];

                     end
                end

            end

            end
    end



        return M_1,M_2,M_3,M_4;

end

