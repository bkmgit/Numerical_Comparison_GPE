function Initial_Value(Simplices,Mesh2Space,Nodes,SpaceSize,N_of_Nodes,N_of_Simplices,Mesh)
    



# --------- INITILIZE MATRICES! -----------------


diag1 = 0; diag2 = 1; diag3 = Mesh-2; diag4 = Mesh-1;
vec1 = ones(SpaceSize)*eps(1/10000);  vec2=vec1[1:SpaceSize-diag2]; vec3 = vec1[1:SpaceSize-diag3]; vec4 = vec1[1:SpaceSize-diag4];
M=spdiagm(0=>vec1)+spdiagm(diag2=>vec2)+spdiagm(diag3=>vec3)+spdiagm(diag4=>vec4)+spdiagm(-diag2=>vec2)+spdiagm(-diag3=>vec3)+spdiagm(-diag4=>vec4) ;
	
# -----------------------------------------------



    Bdary =6.0;
    Boundary = 6.0
    bw =load("Ground_State.jld");
    N = bw["N"];

    A_x1 = -Boundary:12/1201:Boundary	
    A_x2 = -Boundary:12/1201:Boundary	

    A = reshape(bw["U"],N,N)/bw["dx"];
    A = [zeros(N+2)'*0im ; zeros(N)*0im A zeros(N)*0im; zeros(N+2)'*0im];

    itp = interpolate(A, BSpline(Quadratic(Line(OnGrid()))))
    sitp = scale(itp, A_x1, A_x2)
    func(x) = 1.0;

    
    println("ASSEMBLING M")
    M = AssembleMatrix_M(Simplices,Nodes,Mesh2Space, func,M,N_of_Simplices,Quad);
    println("PROJECTING")
    b = L2_Projection(Simplices,Nodes,Mesh2Space,N_of_Simplices,N_of_Nodes,SpaceSize,Quad,sitp,Bdary)

    println("L2-Projection done")
    return M\b;
end



















