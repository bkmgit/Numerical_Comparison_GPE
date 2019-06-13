function Start_Quad(Simplexes,Nodes,Mesh2Space,N_of_Simplexes,U,Quadratic_Poly)

X = [0.0 
     0.5
     1.0
     0.0
     0.5
     0.0];
Y = [0.0 
     0.0
     0.0
     0.5
     0.5
     1.0];
#Help functions
  

#---------------------------

#initilise variables -------

    C_u = zeros(Complex{Float64},3);

#--------------------------------



    for i = 1: N_of_Simplexes


        simplex = Simplexes[:,i];

        idx = Mesh2Space[simplex];
        
        fill!(C_u,0.0);

        for k = 1:3; if(idx[k]!=0) C_u[k] = U[idx[k]]; end  end
        
       #interim = Quadratic_Poly[:,i];

        for p = 1:6
            
            Quadratic_Poly[p,i] = abs( U_h(X[p],Y[p],C_u))^2;
     
        end
       
    end


return Quadratic_Poly
end