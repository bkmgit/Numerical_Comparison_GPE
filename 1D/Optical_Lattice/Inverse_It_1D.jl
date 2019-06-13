
# COMMENTS: CALCULATE EIGENFUCNTION of  (-ϵΔ+V)u + b|u|^2u = λu. N.B! u = U/sqrt(dx). λ is p_k



using SparseArrays
using LinearAlgebra
using JLD



ϵ = 1/2;


function tridiag(a,b,c,n)
    #S[I[k],J[k]] = V[k] sparse(I,J,V[])
    I = a*ones(n-1); J = b*ones(n); K = c*ones(n-1);
    return spdiagm(-1=>I)+spdiagm(0=>J)+spdiagm(1=>K)
end

#
b = 1000.0;

N = 12800;
boundaryX = 16; 
dx = (2*boundaryX)/(N+1); 
#b = β*dx;
β = b/dx  ;  
h_max = 10^4;
v(x) = (1*x^2) + 500*sin(x*pi/4)^2;

x = LinRange(-boundaryX,boundaryX,N+2);

x = x[2:end-1]; 

D_2N = 1/dx^2*tridiag(1,-2,1,N);
I = spdiagm(0=>ones(N))
Δ = D_2N;
D_N = 1/(2*dx)*tridiag(-1,0,1,N);

V = zeros(size(x));
for i = 1:size(x,1)
V[i] = real(v(x[i]));
end
V = spdiagm(0=>V);




Id = spdiagm(0=>ones(N));

A_0 = -ϵ*Δ+V;
Linear = [real(A_0) -imag(A_0); imag(A_0) real(A_0)];


B(R,I) = [spdiagm(0=>R.^2)+spdiagm(0=>I.^2) spzeros(N,N); spzeros(N,N) spdiagm(0=>R.^2)+spdiagm(0=>I.^2)];
A(R,I) = Linear+ β/(R'R+I'I)[1]*B(R,I);
p(R,I) = ([R' I']*A(R,I)*[R;I]/(R'R+I'I)[1])[1];
J(R,I,f) = Linear*f+ β/(R'R+I'I)[1]*( [3*spdiagm(0=>R.^2)+spdiagm(0=>I.^2) 2*spdiagm(0=>R)*spdiagm(0=>I) ; 2*spdiagm(0=>R)*spdiagm(0=>I) spdiagm(0=>R.^2)+3*spdiagm(0=>I).^2]*f)-2/(R'R+I'I)[1]*B(R,I)*[R;I]*dot([R;I],f) ;

u0 = ones(size(x))*10;#exp(-X.^2-Y.^2);
u = (1+2im)*vec(u0);

R = real(u); I = imag(u);
scale = norm([R;I]);
#R = R/scale; I = I/scale;

v_k = [R;I];
v_k = v_k/norm(v_k);

λ = 0;

for i = 1:100;
global R, I, v_k
    p_k = p(R,I);
    
    f_k = p_k*v_k-A(R,I)*v_k;
    g_k = J(R,I,f_k);
    e_k = (-g_k+p_k*f_k)- dot(v_k,(-g_k+p_k*f_k))*v_k + dot(v_k,(A(R,I)*f_k-p_k*f_k))*v_k;
    h_k = (2*ϵ/norm(e_k))^0.5
    if(h_k>h_max) h_k=h_max end
#    if(norm(f_k)<= 10.0^-10) break end;
   
    σ = p_k-1/h_k;

    # DO THE ITERATION!
    #v_k = v_k/norm(v_k);
    println(norm(f_k))
    C = Linear + β*[3*spdiagm(0=>R.^2)+spdiagm(0=>I.^2) 2*spdiagm(0=>R)*spdiagm(0=>I) ; 2*spdiagm(0=>R)*spdiagm(0=>I) spdiagm(0=>R.^2)+3*spdiagm(0=>I).^2]-σ*spdiagm(0=>ones(2*N));
    #C=lufact(C);
    u1 = C\v_k;
    w = 2*β*(C\(B(R,I)*v_k));
    u2 = dot(v_k, u1)/(1-dot(v_k,w))*w;
    
    
    v_k = (u1+u2);
    v_k = v_k/norm(v_k);

    R = v_k[1:N]; I =v_k[N+1:end]

    λ=p_k;
end

println(λ)

U = R+1im*I;
   save("Initial_Value"*".jld","U",U/sqrt(dx))

