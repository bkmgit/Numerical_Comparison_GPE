u(x,t) = (8*exp(4im*t)*(9*exp(-4*x)+16*exp(4*x))-32*exp(16im*t)*(4*exp(-2*x)+9*exp(2*x)))/(-128*cos(12*t)+4*exp(-6*x)+16*exp(6*x)+81*exp(-2*x)+64*exp(2*x))


t = 2.0;

K = 2^11

Mesh = 51200
X = [-20+i*(40/Mesh) for i = 0:Mesh];
h = 40/Mesh;
a = -20.0;
b = 20.0;
mu = 2*pi/(b-a)*(-Mesh/2:Mesh/2-1)


U_ref = [u(x,t) for x in X]
U_ref = U_ref[1:end-1];


b = load("./T2/SP2_Nt_"*string(K)*".jld")
U = b["U"];
U_h = U;


E = U_h-U_ref;

L2 = sqrt(sum(abs.(E).^2*h))

Dx = fftshift(fft(E)).*(1im*mu);

Dx = ifft(ifftshift(Dx));

		
H1 = sqrt(sum(abs.(Dx).^2)*h);

println(L2," ", H1)

