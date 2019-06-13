function AssembleMatrix_J_F_formatFEM_IM( Simplices::Array{Int64,2},Nodes::Array{Float64,2},Mesh2Space::Array{Int64,1} ,U_it::Array{Complex{Float64},1},U0::Array{Complex{Float64},1}, M_1::SparseMatrixCSC{Float64,Int64},M_2::SparseMatrixCSC{Float64,Int64},M_3::SparseMatrixCSC{Float64,Int64},M_4::SparseMatrixCSC{Float64,Int64}, N_of_Simplices::Int64,θ::Float64,Quad::Dict{String,Array{Float64,1}})
#Calculates <|phi|^2 phi_j, phi_i>
    
    
N_gp=( Int(Quad["N_gp"][1]))



    phi = [0.0 0.0 0.0];
    u = 0.0+0.0im;
    u0 = 0.0+0.0im;
    
    C = zeros(3)*0im;
    C_0 = zeros(3)*0im;
   
    loc_1 = zeros(3,3);
    loc_2 = zeros(3,3);
    loc_3 = zeros(3,3);
    loc_4 = zeros(3,3);


    for i = 1: N_of_Simplices
       

        simplex = Simplices[:,i];
        
        idx = Mesh2Space[simplex];

        Vertex = Nodes[:,simplex];
        J = hcat(Vertex[:,2]-Vertex[:,1], Vertex[:,3]-Vertex[:,1]);

        DetJ = abs(det(J));

        max = 0.0; max_0 = 0.0;
        
    fill!(C,0);
    fill!(C_0,0);
    for j = 1:3 if(idx[j]!= 0) C[j] = U_it[idx[j]]; if(max<abs(U_it[idx[j]])) max = abs(U_it[idx[j]]) end   
                                   C_0[j] = U0[idx[j]]; if(max_0<abs(U0[idx[j]])) max_0 = abs(U0[idx[j]]) end end  end

#if(max>10.0^-16 && max_0 >10.0^-16)        

        fill!(loc_1,0);
        fill!(loc_2,0);
        fill!(loc_3,0);
        fill!(loc_4,0);

		for gp = 1:N_gp
		    X = Quad["X"][gp]; Y = Quad["Y"][gp];ω =Quad["ω"][gp];            
		    phi = Φ(X,Y,phi);
		    u = U_h(X,Y,C);
		    u0 = U_h(X,Y,C_0);
		    uθ = θ*u+(1-θ)*u0;
		   
			    for i = 1:3
				for j = 1:3
				    loc_1[i,j] += (2*real(uθ)^2+abs(uθ)^2 )*θ*phi[j]*phi[i]*ω
				    
				    loc_2[i,j] += 2*imag(uθ)*real(uθ)*θ*phi[i]*phi[j]*ω;

				    loc_3[i,j] +=  2*real(uθ)*imag(uθ)*θ*phi[i]*phi[j]*ω

				    loc_4[i,j] +=  (2*imag(uθ)^2+abs(uθ)^2 )*θ*phi[j]*phi[i]*ω;               
				end
			    end
		end
	       
        loc_1 = loc_1*DetJ/2.0;
        loc_2 = loc_2*DetJ/2.0;
        loc_3 = loc_3*DetJ/2.0;
        loc_4 = loc_4*DetJ/2.0;

    
         #Add to Global Matrix



        
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

#end
    end



    return M_1,M_2,M_3,M_4;

end

