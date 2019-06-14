  
BoundaryX = -30.0; BoundaryY=70.0;

Mesh = 25600

bw =load("./Meshes/1dMesh_"*string(Mesh)*".jld")
Nodes = bw["Nodes"];  
Simplices=convert(Array{Int64,2}, bw["Simplices"]);
N_of_Nodes = length(Nodes);
N_of_Simplices = size(Simplices,2);
        
Mesh2Space, SpaceSize = Dirichlet_BC( BoundaryX, BoundaryY, N_of_Nodes, Nodes);


#------------------Calculate Reference solution--------------------
t=10;
u(x) = sqrt(2)*exp(1im*(x/2+3*t/4))*sech(x-t)

U_ref = zeros(Complex{Float64},SpaceSize+2);

idx = 1;

for n in Nodes
global idx
  U_ref[idx] = u(n); idx+=1;
end
U_ref = U_ref[2:end-1];

#------------------------Initialize Matrices---------------------
diag1 = 0; diag2 = 1;
vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2];

M   = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
A   = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);
M_b = spdiagm(0=>vec1)+spdiagm(1=>vec2)+spdiagm(-diag2=>vec2);

#----------------------------------------------------------------


#----------------------Assemble Matrices----------------------
func(x) =1.0;
M   = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);
A   = AssembleMatrix_A(Simplices,Nodes,Mesh2Space,A, N_of_Simplices);
M_b = AssembleMatrix_M_formatFEM( Simplices,Nodes,Mesh2Space,U_ref, M_b, N_of_Simplices,Quad)
   
#----------------------------------------------------------#


println("Name", "	", "L2", "	","	", "L1-Rho", "	","	", "H1","	","	", "CPU [s]");

println("-------------------------------------------------------------");


for Name = ["CN" "IM" "LCN" "RE" "RE2" "TWOSTEP"]

b= load("./T10/"*Name*"_beta-1.0_Nt_256_Mesh_25600.jld")



dec = 10^6;

U = b["U"];

diff = U-U_ref;
L2 = diff'*(M*diff);

rho = abs.(abs.(U).^2-abs.(U_ref).^2);
L2_rho = rho'*(M*rho);
L2_rho = sqrt(L2_rho);

L1_rho = dot(ones(Mesh-1),(M*rho))


L2 = real(sqrt(L2));

H1 = diff'*(A*diff);
H1 = sqrt(real(H1)+L2^2);


println(Name,"	", round(dec*real(L2_rho))/dec, " 	", round(dec*L1_rho)/dec, "	", round(dec*H1)/dec, "	", round(b["Tid"]) ) ;

end



