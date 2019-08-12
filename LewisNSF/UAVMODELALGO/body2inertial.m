function y = body2inertial(u)
vel = u(1:3);
phi = u(4); cp = cos(phi); sp = sin(phi);
the = u(5); ct = cos(the); st = sin(the);
psi = u(6); cs = cos(psi); ss = sin(psi);
mat = [ct*cs sp*st*cs-cp*ss cp*st*cs+sp*ss;...
    ct*ss sp*st*ss+cp*cs cp*st*ss-sp*cs;...
    -st sp*ct cp*ct];
y = mat*vel;