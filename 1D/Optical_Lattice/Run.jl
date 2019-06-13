include("../Schemes/Load_Schemes.jl")

using Interpolations

pathway = "./T1/";

#Constant Declarations

const global β = 1000.0;
const global T = 1.0;
const global ϵ = 1.0;
const global θ = 0.5;
const BoundaryX = -16.0; const BoundaryY=16.0;

v(x) = x^2/4+500*sin(x*pi/4)^2


Mesh = 800
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
Nodes_Initial = LinRange(BoundaryX,BoundaryY,12800+2);
U0 = load("Initial_Value.jld")["U"];

itp_real = interpolate(real([0;U0;0]), BSpline(Cubic(Line(OnGrid()))))
sitp_real = scale(itp_real,Nodes_Initial)
itp_imag = interpolate(imag([0;U0;0]), BSpline(Cubic(Line(OnGrid()))))
sitp_imag = scale(itp_imag,Nodes_Initial)

U0 = [sitp_real(x)+1im*sitp_imag(x) for x in Nodes ];
U0 = U0[2:end-1]

#--------------------------------------------------



N_t = 2^11;
tau = T/N_t;


LCN(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_t,  v, pathway, U0); println("LCN-FEM done")

RE(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_t,  v, pathway, U0); println("RE-FEM done")

TWOSTEP(Simplices, Nodes,  N_of_Simplices, N_of_Nodes, Mesh2Space, SpaceSize, Quad,  tau, N_t,  v, pathway, U0); println("TWOSTEP-FEM done")




