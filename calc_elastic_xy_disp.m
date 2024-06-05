function [dispx, dispy, dispz, ...
          ofx, ofy, ofxe, ofye, ...
          tfx, tfy, tfxe, tfye]      = calc_elastic_xy_disp(sx, sy, ...
                                                            fx1, fy1, fx2, fy2, ...
                                                            D, dip, ...
                                                            ss, ds, ts, Pr, bd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                     %%
%%  calc_elastic_xy_disp.m                             %%
%%                                                     %%
%%  There really isn't much here.  This just converts  %%
%%  fault trace, locking depth, and dip information    %%
%%  into the format used for Okada's calculations,     %%
%%  and then calls the Okada.                          %%
%%                                                     %%
%%  All distances are assumed to be in kilometers.     %%
%%  All dislocation components are assumed to be in    %%
%%  meters.                                            %%
%%                                                     %%
%%  Arguments:                                         %%
%%    sx:          x position of station               %%
%%    sy:          y position of station               %%
%%    fx1:         x position of fault start point     %%
%%    fy1:         y position of fault start point     %%
%%    fx2:         x position of fault end point       %%
%%    fy2:         y position of fault end point       %%
%%    D:           fault locking depth                 %%
%%    dip:         fault dip (in radians)              %%
%%    ss:          strike slip component (km)          %%
%%    ds:          dip slip component (km)             %%
%%    ts:          tensile slip component (km)         %%
%%    Pr:          Poisson's ratio                     %%
%%    bd:          burial depth                        %%
%%                                                     %%
%%  Returned variables:                                %%
%%    dispx:       x displacment at each station       %%
%%    dispy:       y displacment at each station       %%
%%    dispz:       z displacment at each station       %%
%%    ofx1:        x position of fault start point     %%
%%                 bottom corner                       %%
%%    ofy1:        y position of fault start point     %%
%%                 bottom corner                       %%
%%    ofx2:        x position of fault start point     %%
%%                 bottom corner                       %%
%%    ofy2:        y position of fault start point     %%
%%                 bottom corner                       %%
%%                                                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate the paramters for Okada's calculation  %%
%%  Also calculate all of the fault plane corners    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[strike, L, W, ofx, ofy, ofxe, ofye, tfx, tfy, tfxe, tfye] = fault_params_to_okada_form(fx1, fy1, fx2, fy2, dip, D, bd);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate elastic deformation  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  [dispx, dispy, dispz]         = okada_plus_op(ofx, ofy, strike, D, dip, L, W, ss, ds, ts, sx, sy, Pr);
[dispx, dispy, dispz]         = okada_disloc(ofx, ofy, strike, D, dip, L, W, ss, ds, ts, sx, sy, Pr);
