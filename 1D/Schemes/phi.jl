function Φ(X::Float64,phi::Array{Float64}) 
       
    phi[1]=1-X
    phi[2]= X;
    return phi
 end
    