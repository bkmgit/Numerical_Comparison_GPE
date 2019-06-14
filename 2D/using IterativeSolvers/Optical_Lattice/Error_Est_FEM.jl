
include("../Schemes/Load_Schemes.jl")

const global T = 1.0;
const global ϵ = 1.0/2.0;
const global θ = 0.5;
const global beta = 2300.0
const BoundaryX = 6.0; const BoundaryY = 6.0;

Mesh = 400;
SpaceSize = 399*399

func(x) = 1.0;

diag1 = 0; diag2 = 1; diag3 = Mesh-2; diag4 = Mesh-1;
vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2]; vec3 = vec1[1:SpaceSize-diag3]; vec4 = vec1[1:SpaceSize-diag4];

M=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;
    
A = spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(-diag2=>vec2)+spdiagm(-diag4=>vec4)+spdiagm(diag4=>vec4);



pathway = "./T1/"


Ref_Sol = load(pathway*"RE_beta_2300.0_Nt_64_Mesh_400.jld")
Ref_Mesh = load("./Meshes/Mesh_400x400.jld")


U_ref = Ref_Sol["U"]; 

U_ref = reshape(U_ref,SpaceSize,1);
U_ref = U_ref[:]



Mesh2Space=Ref_Sol["Mesh2Space"];
Nodes = Ref_Mesh["Nodes"]; Simplices = Ref_Mesh["Simplices"];

Simplices=convert(Array{Int64,2}, Ref_Mesh["Simplices"]);
N_of_Nodes = convert(Int64,Ref_Mesh["N_of_Nodes"]);
N_of_Simplices=convert(Int64,Ref_Mesh["N_of_Simplices"]);

M = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);
A = AssembleMatrix_A(Simplices,Nodes,Mesh2Space,A, N_of_Simplices);

Decimals = 10.0^6

Sol = load(pathway*"LCN_beta_2300.0_Nt_64_Mesh_400.jld")

U = Sol["U"];

diff = U-U_ref;
rho = abs.(abs.(U).^2-abs.(U_ref).^2)
L1_rho = dot(ones(SpaceSize),M*rho);
L2_rho = sqrt(dot(rho,M*rho));
L2 = sqrt(abs(diff'*M*diff));
H1 = sqrt(abs(diff'*A*diff));


println("L2: ",round(Decimals*(L2))/Decimals, "	H1: ", round(Decimals*(H1))/Decimals,"	L1_rho: ", round(Decimals*(L1_rho))/Decimals);


