using Plots
using JLD

t = 200.0;
plot()


u(x,t) =  (8*exp(4im*t)*(9*exp(-4*x)+16*exp(4*x))-32*exp(16im*t)*(4*exp(-2*x)+9*exp(2*x)))/(-128*cos(12*t)+4*exp(-6*x)+16*exp(6*x)+81*exp(-2*x)+64*exp(2*x)) ;


X = LinRange(-20,20,1*51200+1);
X = X[2:end-1];

for k = 2 .^[5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]

#	b = load("./T200_Finer/Relaxation_beta-2.0_tau"*string(t/k)* "_Mesh102400.jld");
# b = load("./T200_Finer/Crank_beta-2.0_tau"*string(t/k)*"_Mesh102400.jld");

 b = load("./T200/Crank_beta-2.0_tau"*string(t/k)*"_Mesh51200.jld");
	U = b["U"];
	Energy = b["Energy"];

	#plot(abs.(U).^2);
	#savefig("./T500/Re_FEM_"*string(k)*".png");

	#plot()
	CPU = 6.33*k/2^17;
	CPU = round(100*CPU)/100
#	plot([0; Energy[1:end-1,2]],[-48.0; Energy[1:end-1,1]],label = string(k)*" CPU Time ~"*string(CPU));

	plot(X,abs.(U).^2,label = string(k)*" CPU Time ~"*string(CPU));
	plot!(X, [abs(u(x,200)).^2 for x in X],label="ref")
	savefig("Abs_FEMCN_"*string(k)*".png");
end

	#savefig("./T500/Re_FEM.png");



#=
x = LinRange(-5,5,1000);
idx = 1
for t = 0:0.01:1;
global idx
	U = [abs(u(X,t))^2 for X in x];
	plot(x,U);savefig("./Phases/abs_"*string(idx)*".png"); idx+=1
end
=#
