function Φ(X::Float64,Y::Float64,phi::Array{Float64}) 
       
    phi[1]=1-X-Y
    phi[2]= X;
    phi[3]= Y;
    return phi
 end
    
U_h(X::Float64,Y::Float64,C::Array{Complex{Float64},1}) = (1-X-Y)*C[1] +C[2]*X +C[3]*Y ;
