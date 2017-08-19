function [ux, uz] = DipTwoD(xp1, zp1, xp2, zp2, xobs, ds, ts)
% Convert xp1, yp1, xp2, and yp2 parameters into Okada style parameters
delta                         = atan2(zp2-zp1, xp1-xp2);
fx1                           = xp1 + zp1/tan(delta);
fx2                           = fx1;
D                             = zp2;
bd                            = zp1;

% Fixed values
fy1                           = 5000;
fy2                           = -5000;

%  Calculate surface displacements
[ux, ~, uz, ~, ~, ~, ~, ~, ~, ~, ~]    = calc_elastic_xy_disp(xobs, zeros(size(xobs)), fx1, fy1, fx2, fy2, D, delta, 0, ds, ts, 0.25, bd);

