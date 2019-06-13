using FFTW

function SP2_FFT(a,b,N,Nodes,beta,epsilon,U,V_func,T,K,pathway)
 Tid = time();   
    #-------------------------------------------------------------#

    Nodes = Nodes[1:end-1];
    
    V = [V_func(n) for n in Nodes];
    
    #--------------------------------------------------------------#
    
    k = T/K;
  #  print(k)
    
    mu = 2*pi/(b-a)*(-N/2:N/2-1)

#---------------------------ENERGY-Measurements----------------------------#
N_E = 200;
if(K>N_E)
	E =[1+K/N_E*i for i =0:N]
	E = round.(E)
	Energy = zeros(N_E+1,2);
else
	E = [1+i for i = 0:K]
	E = round.(E)
	Energy = zeros(Int(K)+1,2);
end 

idx = 1
#--------------------------------------------------------------------------#

it = 1
    #----------------------TIME STEPPING ---------------------------#
    for i = 1:K
       # Evolve without laplacian (HALF STEP)
	if( i == E[idx]) 
		Dx = fftshift(fft(U)).*(1im*mu);

		Dx = ifft(ifftshift(Dx));

		dx = 40/N;
		energy = epsilon*(abs.(Dx).^2)+ V.*(abs.(U).^2) + beta/2*abs.(U).^4
		
		Energy[idx,1] = sum(energy)*dx;
		Energy[idx,2] = (i-1)*k;
		
		idx += 1
	end
	Ustar = exp.(-1im*(V+beta*abs.(U).^2)*k/2).*U
	#println(i*k)


        Ustarhat = fftshift(fft(Ustar));
	
	#Laplace and time
	Ustarhat = Ustarhat.*exp.(-1im*epsilon*2*mu.^2*k/2)
	U2star = ifft(ifftshift(Ustarhat))


        #Evolve with laplace

	

        #Evolve without laplacian (HALF STEP)
         U = exp.(-1im*(V+beta*abs.(U2star).^2)*k/2).*U2star 
	
    end
    #------------------------End--------------------------------#
		Dx = fftshift(fft(U)).*(1im*mu);

		Dx = ifft(ifftshift(Dx));

		dx = 40/N;
		energy = epsilon*(abs.(Dx).^2)+ V.*(abs.(U).^2) + beta/2*abs.(U).^4
		
		Energy[end,1] = sum(energy)*dx;
		Energy[end,2] = T;    


    Tid = time()-Tid;
   save(pathway*"SP2_Nt_"*string(K)*".jld","U",U,"Energy",Energy,"Tid",Tid)
	return U
   
end




