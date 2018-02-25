function [ K1,K2 ] = PLL_Design( Zeta, deltaF, Ts )

N = 1;

Bn = deltaF / (2*pi*sqrt(2)*Zeta)

K1 = (4*Zeta/N)*(Bn*Ts/(Zeta+(1/(4*Zeta))));
K1 = K1 / (1+(2*Zeta/N)*(Bn*Ts/(1/(4*Zeta)))+(Bn*Ts/(N*(Zeta+(1/(4*Zeta)))))^2)

K2 = (4/(N^2)) * ((Bn*Ts/(Zeta+(1/(4*Zeta))))^2);
K2 = K2 / (1+(2*Zeta/N)*(Bn*Ts/(1/(4*Zeta)))+(Bn*Ts/(N*(Zeta+(1/(4*Zeta)))))^2)

end

