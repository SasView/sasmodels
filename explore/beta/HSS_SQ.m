function Sq=HSS_SQ(volF,Q)
%Sq=HSS_SQ(volF,Q), Q is a normalzied scattering wavevector by the diamter of the particle

a1=-(1+2*volF)^2/(1-volF)^4;
a2=6*volF*(1+0.5*volF)^2/(1-volF)^4;
a3=volF*a1/2;
%pi = 3.14159265358979;

cq=4*pi*(a1./Q.^3 .* (sin(Q)-Q.*cos(Q)) + ...
    a2./Q.^4 .* (2*Q.*sin(Q)- (Q.^2 - 2) .* cos(Q) - 2) + ...
    a3./Q.^6 .* ((4*Q.^3 - 24 * Q) .* sin(Q) - (Q.^4 -12 * Q.^2 + 24) .* cos(Q) + 24) );

Sq = 1 ./ (1 - cq * volF * 6 / pi);