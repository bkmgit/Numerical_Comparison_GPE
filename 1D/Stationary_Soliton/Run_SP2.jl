using JLD

include("../SP2/SP2_FFT.jl")


pathway = "./T2/";

const global β = -2.0;
const global T = 2.0;
const global ϵ = 1.;
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


# ----------------------------Calculate U0 -----------------------------------
t= 0.0;
u0(x) = (8*exp(4im*t)*(9*exp(-4*x)+16*exp(4*x))-32*exp(16im*t)*(4*exp(-2*x)+9*exp(2*x)))./(-128*cos(12*t)+4*exp(-6*x)+16*exp(6*x)+81*exp(-2*x)+64*exp(2*x))

U0 = zeros(Complex{Float64},Mesh);
for i = 1:Mesh
	U0[i] = u0(Nodes[i]) 
end
#--------------------------------------------------------------------------------           
 

N_t = 2^11
U = SP2_FFT(BoundaryX,BoundaryY,Mesh,Nodes,β,ϵ,U0,v,T,N_t,pathway); println("SP2 done");




















