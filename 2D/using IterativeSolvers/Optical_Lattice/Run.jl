include("Initial_Value.jl")
include("../Schemes/Load_Schemes.jl")

pathway = "./T1/";


const global β = 2300.0;
const global T = 1.;
const global ϵ = 1.0/2;
const global θ = 0.5;
const BoundaryX = -6.0; const BoundaryY=6.0;
const global Omega = 0.0

gammax = 4;
gammay = 8;

v(x)  = ((gammax*x[1])^2+(gammay*x[2])^2)/2 +1000*( (abs(x[1])-4.5)^5 *(abs(x[1])>=4.5) +(abs(x[2])-4.5)^5 *(abs(x[2])>=4.5) )+787*(sin(pi*x[1]/2)^2 +sin(pi*x[2]/2).^2);



	Mesh = 400;
	bw = load("./Meshes/Mesh_400x400.jld")	
	Nodes = bw["Nodes"];  
	Simplices=convert(Array{Int64,2}, bw["Simplices"]);
	N_of_Nodes = convert(Int64,bw["N_of_Nodes"]);
	N_of_Simplices=convert(Int64,bw["N_of_Simplices"]);
        #---------------Space smaller than mesh--------------


	Mesh2Space, SpaceSize = Define_Space( BoundaryX, BoundaryY, N_of_Nodes,Nodes)
	
#	U0 = Initial_Value(Simplices,Mesh2Space,Nodes,SpaceSize,N_of_Nodes,N_of_Simplices,Mesh)
	
b = load("./Ground_State/Ground_State.jld")
U0 = b["U"] 
   
N_t = 2^6
tau = T/N_t

    
RE(Simplices,Nodes,N_of_Simplices,N_of_Nodes,Mesh2Space,SpaceSize,tau,N_t,v,pathway,U0); println("RE-FEM done")
LCN(Simplices,Nodes,N_of_Simplices,N_of_Nodes,Mesh2Space,SpaceSize,tau,N_t,v,pathway,U0); println("LCN-FEM done")
#CN(Simplices,Nodes,N_of_Simplices,N_of_Nodes,Mesh2Space,SpaceSize,tau,N_t,v,pathway,U0)
#IM(Simplices,Nodes,N_of_Simplices,N_of_Nodes,Mesh2Space,SpaceSize,tau,N_t,v,pathway,U0)



