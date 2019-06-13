function tridiag(a,b,c,n)
    #S[I[k],J[k]] = V[k] sparse(I,J,V[])
    I = a*ones(n-1); J = b*ones(n); K = c*ones(n-1);
    return spdiagm(-1=>I,0=>J,1=>K )
    
end


function linspace(a,b,N)
	return  [i*(b-a)/(N-1)+a for i = 0:N-1]
end







function Doit(v,b,E_0)
Ω = 0*.8;
ϵ = 1.0/2;
N = 399;  #NUMBER OF INNER POINTS!
boundaryX = 6; boundaryY = 6;
dx = (2*boundaryX)/(N+1); dy = (2*boundaryY)/(N+1);
β = b/dx^2;
#b = β*dx^2;
h_max = 10^0;


x = linspace(-boundaryX,boundaryX,N+2);
y = linspace(-boundaryY,boundaryY,N+2);

x = x[2:end-1]; y = y[2:end-1];

D_2N = 1/dx^2*tridiag(1,-2,1,N);
#I = eye(N);
I = spdiagm(0=>ones(N));#Matrix(I,N,N);
Δ = kron(D_2N,I) + kron(I,D_2N);
D_N = 1/(2*dx)*tridiag(-1,0,1,N);
rot = kron(spdiagm(0=>y),D_N)-kron(D_N,spdiagm(0=>x));

X = kron(x,ones(1,N));
Y = kron(y',ones(N,1));
V = v(X,Y);




Id = spdiagm(0=>ones(N*N));

A_0 = -ϵ*Δ+1im*Ω*rot+ spdiagm(0=>vec(V));
Linear = [real(A_0) -imag(A_0); imag(A_0) real(A_0)];


B(R,I) = [spdiagm(0=>R.^2)+spdiagm(0=>I.^2) spzeros(N*N,N*N); spzeros(N*N,N*N) spdiagm(0=>R.^2)+spdiagm(0=>I.^2)];
A(R,I) = Linear+ β/(R'R+I'I)[1]*B(R,I);
p(R,I) = ([R' I']*A(R,I)*[R;I]/(R'R+I'I)[1])[1];
J(R,I,f) = Linear*f+ β/(R'R+I'I)[1]*( [3*spdiagm(0=>R.^2)+spdiagm(0=>I.^2) 2*spdiagm(0=>R)*spdiagm(0=>I) ; 2*spdiagm(0=>R)*spdiagm(0=>I) spdiagm(0=>R.^2)+3*spdiagm(0=>I).^2]*f)-2/(R'R+I'I)[1]*B(R,I)*[R;I]*dot([R;I],f) ;



    Omega = 0.8;
    phi0(X,Y) = 1/sqrt(pi)*exp.(-(X.^2+Y.^2)./2)
    u0(X,Y) = (1-Omega)*phi0(X,Y).+Omega*X.*phi0(X,Y).+1im*Omega*Y.*phi0(X,Y); #+exp(-2*((x-2)^2+y^2))+exp(-2*((x+2)^2+(y+2)^2));




#Bound = ((X.^2+Y.^2).<20);
#ue=ue.*Bound;
u = vec(u0(X,Y));

#u += 2*rand(N,N)+1im*rand(N,N);
#u = vec(u);
R = real(u); I = imag(u);
scale = norm([R;I]);
#R = R/scale; I = I/scale;


v_k = [R;I];
v_k = v_k/norm(v_k);


for h_max = [10000 1000 1000 500 500 500 ]
for i = 1:10;
    p_k = p(R,I);
    
    f_k = p_k*v_k-A(R,I)*v_k;
    g_k = J(R,I,f_k);
    e_k = (-g_k+p_k*f_k)- dot(v_k,(-g_k+p_k*f_k))*v_k + dot(v_k,(A(R,I)*f_k-p_k*f_k))*v_k;
    h_k = (2*ϵ/norm(e_k))^0.5
    if(h_k>h_max) h_k=h_max end
   
    σ = p_k-1/h_k;

    # DO THE ITERATION!
    #v_k = v_k/norm(v_k);
println(round(norm(f_k)*10^7)/10^7," ", p_k," ", h_max)
   
    if(norm(f_k)<10.0.^-12) break end

    C = Linear + β*[3*spdiagm(0=>R.^2)+spdiagm(0=>I.^2) 2*spdiagm(0=>R)*spdiagm(0=>I) ; 2*spdiagm(0=>R)*spdiagm(0=>I) spdiagm(0=>R.^2)+3*spdiagm(0=>I).^2]-σ*spdiagm(0=>ones(2*N^2));
    C=lu(C);
    u1 = C\v_k;
    w = 2*β*(C\(B(R,I)*v_k));
    u2 = dot(v_k, u1)/(1-dot(v_k,w))*w;
    
    
    v_k = (u1+u2);
    v_k = v_k/norm(v_k);

    R = v_k[1:N^2]; I =v_k[N^2+1:end]

end

end

U = R+1im*I;
#file = matopen(join(["Ground_State_Vh_2_p_25_E",string(E_0),"_b",string(b),".mat"]),"w");
#write(file,"U",U)
#write(file,"dx",dx)
#write(file,"b",b)
#write(file,"N",N)
#close(file)

save("Ground_State_invIt.jld","U",U,"dx",dx,"b",b,"N",N)
end
















# p_0  = 2.577

p_0 = 2.5;

# A_0 = 10^5 * 1.05 * pi / ( 8.53 * 8.53 * 1.441 * 24.0);

 A_0 = 100

for E_0 = [ 10]
#E_0 = 11.0;


v(x,y)  = (x.^2+y.^2)/2 .+1000*( (abs.(x).-4.5)^5 .*(abs.(x).>=4.5) .+(abs.(y).-4.5).^5 .*(abs.(y).>=4.5) ).+787*(sin.(pi*x/2).^2 .+sin.(pi*y/2).^2);

#b = 200.0;
b = 2300#1163.115924633;  # NON-LINEARITY IN EQUATION!

Doit(v,b,E_0)

end





