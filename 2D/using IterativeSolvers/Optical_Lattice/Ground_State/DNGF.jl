#using Interpolations
using IterativeSolvers
using SparseArrays
using LinearAlgebra
using JLD

include("Initial_Value.jl")
include("../../Schemes/Define_Space.jl")
#include("E_Norm.jl")

include("../../Schemes/AssembleMatrix_M.jl")
include("../../Schemes/AssembleMatrix_M_formatFEM.jl")
include("../../Schemes/AssembleMatrix_A.jl")
include("../../Schemes/AssembleMatrix_Lz.jl")
include("../../Schemes/phi.jl")
include("../../Schemes/Quadrature_7.jl")
include("../../Schemes/L2_Projection.jl")

include("../../Schemes/Reinit_M.jl")

BoundaryX = 6.0;BoundaryY=6.0;

gammax = 1.0;
gammay = 1.0
v(x)  = ((gammax*x[1])^2+(gammay*x[2])^2)/2 +1000*( (abs(x[1])-4.5)^5 *(abs(x[1])>=4.5) +(abs(x[2])-4.5)^5 *(abs(x[2])>=4.5) )+787*(sin(pi*x[1]/2)^2 +sin(pi*x[2]/2).^2);


Mesh = 400;

bw = load(join(["../Meshes/", "Mesh_",string(Mesh),"x",string(Mesh),".jld"]))
Nodes = bw["Nodes"];  
Simplices=convert(Array{Int64,2}, bw["Simplices"]);
N_of_Nodes = convert(Int64,bw["N_of_Nodes"]);
N_of_Simplices=convert(Int64,bw["N_of_Simplices"]);
        #---------------Space smaller than mesh--------------


Mesh2Space, SpaceSize = Define_Space( BoundaryX, BoundaryY, N_of_Nodes,Nodes)
bw = load("Ground_State.jld")
U = bw["U"]








function DNFG(Simplices,Nodes,Mesh2Space,SpaceSize,N_of_Simplices,N_of_Nodes,U,Mesh)





    #INITIALIZE SPARSE MATRICES
    diag1 = 0; diag2 = 1; diag3 = Mesh-2; diag4 = Mesh-1;
    vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2]; vec3 = vec1[1:SpaceSize-diag3]; vec4 = vec1[1:SpaceSize-diag4];


    M=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

    A=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(-diag2=>vec2)+spdiagm(-diag4=>vec4)+spdiagm(diag4=>vec4);


    M_β=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

        M_v=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

    M_L=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;

    func(x)=1.0;
    #BUILD MATRICES
    M = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);
    A=AssembleMatrix_A(Simplices,Nodes,Mesh2Space,A, N_of_Simplices);
    M_v = AssembleMatrix_M(Simplices,Nodes,Mesh2Space,v,M_v, N_of_Simplices,Quad);
    M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  U[:],M_β, N_of_Simplices,Quad ); 
   

beta = 2300.0;
    tau = 1.1;
   Id = spdiagm(0=>ones(SpaceSize));
u_r = zeros(Complex{Float64},SpaceSize);
    	Mass = sqrt(real(U'*(M*U)));

        U = U./(Mass[1])
    for n = 1:100
	 U_old = U;	
	
	Matris = Id + tau*(0.5*A+M_v+beta*M_β);
        U = U+tau*(-0.5*A*U-M_v*U-beta*M_β*U);
	#U = Matris\U;

	RHS = U;
	println("SOLVING GMRES")
        #tic()
        r = RHS-Matris*U;
        bicgstabl!(u_r, Matris, r[:],tol=10.0^-10);
        U+=u_r;        

	Mass = sqrt(real(U'*(M*U)));
        U = U./(Mass[1])
	
       M_β = Reinit_M(Simplices,Nodes,Mesh2Space, M_β, N_of_Simplices);
       M_β = AssembleMatrix_M_formatFEM(Simplices,Nodes,Mesh2Space,  U[:],M_β, N_of_Simplices,Quad ); 
          println(n)  


#
	Energy = (U'*(0.5*A+M_v+beta*M_β/2)*U)
	Chem_Pot = Energy + U'*beta/2*(M_β*U)
	println("CHEM_POT"," ",real(Chem_Pot))

	println("ENERGY"," ", round(real(Energy)*10^7)/10^7);
	println("norm(Diff):" ," ", norm(U-U_old));


    end
save("Ground_State.jld","U",U)
end



DNFG(Simplices,Nodes,Mesh2Space,SpaceSize,N_of_Simplices,N_of_Nodes,U,Mesh)
