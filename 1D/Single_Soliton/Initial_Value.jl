
function Initial_Value(Mesh2Space::Array{Int64,1},Nodes::Array{Float64,1},SpaceSize::Int64,N::Int64)        

  α=1.0; q = 1.0; c=1.0; 
 u0(x) =√(2*α/q)*exp(0.5*1im*c*x)*sech(x*√(α))


            index = 1;
            U0 = zeros(SpaceSize)*0im;
            for i = 1:N
                if(Mesh2Space[i]!= 0)
                    x = Nodes[i]; 
                    U0[index] = u0(x);
                    index = index+1;
                end
            end
         return   U0 = U0[:];
end




