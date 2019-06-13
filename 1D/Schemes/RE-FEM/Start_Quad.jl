function Start_Quad(Simplexes,Nodes,Mesh2Space,N_of_Simplexes,U,Quadratic_Poly)

X = [0.0 
     1.0
     0.5];

#Help functions
  

#---------------------------

#initilise variables -------

    C_u = zeros(Complex{Float64},2);

#--------------------------------
U_h(X,C) = (1-X)*C[1]+X*C[2];


    for i = 1: N_of_Simplexes


        simplex = Simplexes[:,i];

        idx = Mesh2Space[simplex];
        
        fill!(C_u,0.0);

        for k = 1:2; if(idx[k]!=0) C_u[k] = U[idx[k]]; end  end
        
       #interim = Quadratic_Poly[:,i];

        for p = 1:3
            
            Quadratic_Poly[p,i] = abs( U_h(X[p],C_u))^2;
     
        end
       
    end


return Quadratic_Poly
end
