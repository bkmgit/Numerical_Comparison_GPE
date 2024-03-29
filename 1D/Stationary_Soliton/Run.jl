include("../Schemes/Load_Schemes.jl")


include("Initial_Value.jl");


pathway = "./T2/";

#Constant Declarations

const global β = -2.0;
const global T = 2.0;
const global ϵ = 1.0;
const global θ = 0.5;
const BoundaryX = -20.0; const BoundaryY=20.0;

v(x)  = 0;


Mesh = 51200
#--------------- Load Mesh--------------- 
bw =load("./Meshes/1dMesh_"*string(Mesh)*".jld")
Nodes = bw["Nodes"];  
Simplices=convert(Array{Int64,2}, bw["Simplices"]);
N_of_Nodes = length(Nodes);
N_of_Simplices = size(Simplices,2);
#-----------------------------------------

#---------------Dirichlet BC --------------
Mesh2Space, SpaceSize = Dirichlet_BC( BoundaryX, BoundaryY, N_of_Nodes,Nodes);
#-----------------------------------------------------------


#-------------Interpolation operator--------------- 
U0= Initial_Value(Mesh2Space,Nodes,SpaceSize,N_of_Nodes) 
#--------------------------------------------------


N_t = 2^7;

tau = T/N_t;


LCN(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_t,  v, pathway, U0); println("LCN-FEM done")

RE(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_t,  v, pathway, U0); println("RE-FEM done")

RE2(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_t,  v, pathway, U0); println("RE2-FEM done")

CN(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_t,  v, pathway, U0); println("CN-FEM done")

IM(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_t,  v, pathway, U0); println("IM-FEM done")

TWOSTEP(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_t,  v, pathway, U0); println("TWOSTEP-FEM done")














