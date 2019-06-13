x_low = -6; x_high = 6;
y_low = -6; y_high = 6;

function linspace(a,b,N)
   return [a+i*(b-a)/(N-1) for i = 0:N-1] 
end

for N = 500
Nx = N;  Ny = N;

x = linspace(x_low,x_high,Nx+1);
y = linspace(y_low,y_high,Ny+1);


Nodes = zeros(2,(Nx+1)*(Ny+1));
Simplexes = zeros(3,2*Nx*Ny);

for j = 1:Ny

    for i = 1:Nx
      Nodes[:,i+(j-1)*(Nx+1)] = [x[i];y[j]]; 
        
      Simplexes[:,2*(i+(j-1)*Nx)-1] = [i+(j-1)*(Nx+1); i+1+(j-1)*(Nx+1); i+j*(Nx+1) ];
      Simplexes[:,2*(i+(j-1)*Nx)] = [i+1+(j-1)*(Nx+1); i+1+(j)*(Nx+1); i+j*(Nx+1) ];
        
    end
    i = Nx;
    Nodes[:,i+(j-1)*(Nx+1)+1]=[x_high,y[j]];
end

j  = Ny+1;
    for i = 1:Nx
        Nodes[:,i+(j-1)*(Nx+1)] = [x[i];y[j]]; 
    end
    i = Nx;
    Nodes[:,i+(j-1)*(Nx+1)+1]=[x_high,y[j]];
    
N_of_Nodes = (Nx+1)*(Ny+1);
N_of_Simplexes = 2*Nx*Ny;

save("Mesh_"*string(Nx)*"x"*string(Ny)*".jld", "Nodes", Nodes, "Simplexes",Simplexes, "N_of_Nodes",N_of_Nodes,"N_of_Simplexes",N_of_Simplexes)
#file = matopen("matfile.mat")    
#save(join(['Mesh_',num2str(Nx),'x',num2str(Ny)]),'Nodes','Simplexes','N_of_Nodes','N_of_Simplexes')
end
