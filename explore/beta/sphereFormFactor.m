function data=sphereFormFactor(Q,radius)

data= ( 3 * (sin(Q*radius) - Q*radius.*cos(Q*radius)) ./ (Q*radius).^3).^2 ;

return;