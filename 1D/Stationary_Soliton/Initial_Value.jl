function Initial_Value(Mesh2Space::Array{Int64,1},Nodes::Array{Float64,1},SpaceSize::Int64,N::Int64)        
#=
            α=0.5; q = 1.0; c₁=1.0; c2 = 0.1; δ=25;
    
            u0(x) = √(2*α/q)*exp(0.5*1im*c₁*x)*sech(x*√(α))+exp(0.5*1im*c2*(x-δ))*sech((x-δ)/√(α)) Two solitons
    =#
 #   α=10.0; q = 1.0; c=1.0; 
   # u0(x) =√(2*α/q)*exp(0.5*1im*c*x)*sech(x*√(α))
    t= 0.0;
  u0(x) = (8*exp(4im*t)*(9*exp(-4*x)+16*exp(4*x))-32*exp(16im*t)*(4*exp(-2*x)+9*exp(2*x)))./(-128*cos(12*t)+4*exp(-6*x)+16*exp(6*x)+81*exp(-2*x)+64*exp(2*x))
   #a = 10.0;

            index = 1;
            U0 = zeros(SpaceSize)*0im;
            for i = 1:N
                if(Mesh2Space[i]!= 0)
                    x = Nodes[i]; 
                    #	if(abs(x)==8)
                    #		U0[index] =0.0;
                    #	elseif(abs(y)==8)
                    #		U0[index]=0.0;
                    #	else	
                    U0[index] = u0(x);
                    #	end
                    index = index+1;
                end
            end
         return   U0 = U0[:];
end
