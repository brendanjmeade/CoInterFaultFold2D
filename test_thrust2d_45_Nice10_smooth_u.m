close all; clear all;

% To do:
% - Graphically represent slip on fault and fold
% - Interseismic
% - Unify fonts (to override startup.m)

% Declare station positions
n_stations                    = 1000;
xmin                          = -30;
xmax                          = 30;
sx                            = linspace(xmin, xmax, n_stations);
nlay                          = 7;
FontName                      = 'Helvectica';

% Primary fault surface
faultx1                       = 7;
faultz1                       = 0;
faultx2                       = 0;
faultz2                       = 5;
faultm                        = -(faultz2-faultz1)/(faultx2-faultx1);
gamma                         = abs(atand(-(faultz2-faultz1)/(faultx2-faultx1)));
R                             = [cosd(gamma) -sind(gamma) ; sind(gamma) cosd(gamma)]; % rotation matrix for 
bisectorAngle                 = (180-gamma)/2;
[uxfault, uzfault]            = DipTwoD(faultx1, faultz1, faultx2, faultz2, sx, 1, 0);

% Primary fold surface
foldx1                        = 0;
foldz1                        = 5;
% foldx2                        = -3.5; % Just an arbitrary example
foldx2                        = -abs(faultz2-faultz1)/tand(bisectorAngle); % With bisector
foldz2                        = 0;
foldm                         = -(foldz2-foldz1)/(foldx2-foldx1);
x1                            = faultx2-faultx1;
z1                            = faultz2-faultz1;
x2                            = foldx2-foldx1;
z2                            = foldz2-foldz1;
dsfold                        = x1*x2 + z1*z2; % dip-slip component on fold
tsfold                        = -(1-dsfold); % tensile component on fold
magfold                       = abs(dsfold)+abs(tsfold);
dsfold                        = dsfold/magfold;
tsfold                        = tsfold/magfold;
[uxfold_ds, uzfold_ds]        = DipTwoD(foldx1, foldz1, foldx2, foldz2, sx, dsfold, 0);
[uxfold_ts, uzfold_ts]        = DipTwoD(foldx1, foldz1, foldx2, foldz2, sx, 0, tsfold);

% Calls to classical Okada for the horizontal case only.
xdetach                       = 8.6023; % Same length as ramp...at least for now
                                        %xdetach                       = 0; % Same length as ramp...at least for now

[uxdetach, ~, uzdetach]       = okada_plus_op(0, -5000, pi/2, 5, 0, 10000, xdetach, 0, -1, 0, sx, zeros(size(sx)), 0.25);

% Unit block displacment
idx1                          = find(sx < faultx1);
idx2                          = find(sx > foldx2);
idx                           = intersect(idx1, idx2);
temp                          = R*[1 ; 0];
uxblock                       = zeros(size(sx));
uzblock                       = zeros(size(sx));
uxblock(idx)                  = temp(1);
uzblock(idx)                  = temp(2);

% Total displacements
uxfault                       = uxfault(:);
uxfold_ds                     = uxfold_ds(:);
uxfold_ts                     = uxfold_ts(:);
uxdetach                      = uxdetach(:);
uxblock                       = uxblock(:);
uzfault                       = uzfault(:);
uzfold_ds                     = uzfold_ds(:);
uzfold_ts                     = uzfold_ts(:);
uzdetach                      = uzdetach(:);
uzblock                       = uzblock(:);

UX                            = uxfault + uxfold_ds + uxfold_ts + uxdetach + uxblock;
UZ                            = uzfault + uzfold_ds + uzfold_ts + uzdetach + uzblock;

% Get the long term displacements
uxlong                        = uxblock;
uxlong(sx < foldx2)           = 1;
uzlong                        = uzblock;


% Plot displacements (block)
    UX_temp = uxblock;
    UZ_temp = uzblock;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('block');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
    %fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
    %set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
%     plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);

    
% Plot displacements (coseismic fault)
    UX_temp = uxfault;
    UZ_temp = uzfault;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('ramp');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
    plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
%     plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);

    
% Plot displacements (fold tensile)
    UX_temp = uxfold_ts;
    UZ_temp = uzfold_ts;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('fold(tensile)');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
    plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    % plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);

    
    
% Plot displacements (fold tensile)
    UX_temp = uxfold_ds;
    UZ_temp = uzfold_ds;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('fold(dip)');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
    plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    % plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);

% Plot displacements (detachment)
    UX_temp = uxdetach;
    UZ_temp = uzdetach;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('detachment');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);
    
% Plot displacements (ramp + block)
    UX_temp = uxfault + uxblock;
    UZ_temp = uzfault + uzblock;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('ramp + block');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
    %fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
    %set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
    plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    % plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);

    
    
    
    
    
    
    
    
% Plot displacements (fold + block)
    UX_temp = uxfold_ts + uxfold_ds + uxblock;
    UZ_temp = uzfold_ts + uzfold_ds + uzblock;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('fold + block');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
    %fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
    %set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
    plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    % plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
% Plot displacements (fault + detachment)
    UX_temp = uxfault + uxdetach;
    UZ_temp = uzfault + uzdetach;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('ramp + detachment');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
    plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);
    
    
    
    
% Plot displacements (fault + fold)
    UX_temp = uxfault + uxfold_ds + uxfold_ts;
    UZ_temp = uzfault + uzfold_ds + uzfold_ts;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('ramp + fold(tensile) + fold(dip)');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
    plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
    plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    % plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);

    
% Plot displacements (fault + fold + detachment)
    UX_temp = uxfault + uxfold_ds + uxfold_ts + uxdetach;
    UZ_temp = uzfault + uzfold_ds + uzfold_ts + uzdetach;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('ramp + fold(tensile) + fold(dip) + detachment');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
    plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
    plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);
    
% Plot displacements (fault + fold + block)
    UX_temp = uxfault + uxfold_ds + uxfold_ts + uxblock + uxdetach;
    UZ_temp = uzfault + uzfold_ds + uzfold_ts + uzblock + uzdetach;
    UX_total_all = UX_temp;
    UZ_total_all = UZ_temp;
    
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('ramp + fold(tensile) + fold(dip) + block + detachment');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, 'u_z', 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
    %fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
    %set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
    plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
    plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);

    
    
% Plot displacements (fault + fold + block vs. scaled fault : the aliasing problem)
    UX_temp = uxfault + uxfold_ds + uxfold_ts + uxblock + uxdetach;
    UX_temp = uxfault + 1.26.*uxfold_ds + 0.0.*uxfold_ts + 0.0.*uxblock + uxdetach;

    UZ_temp = 1.0*(uxfault+uxdetach);
    scale = UZ_temp(:) \ UX_temp(:);
    UZ_temp = scale*UZ_temp;
    
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('horizontal displacement aliasing');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-w', 'Linewidth', 3);
    plot(sx, UZ_temp, '--b', 'Linewidth', 3);    
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-w', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '--b', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x (all)', 'FontSize', 14);
    text(-22, -1.3-0.03, sprintf('u_x %2.1f(ramp+detachment)', scale), 'FontSize', 14);

    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
    %fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
    %set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
    plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
    plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);

    
% Plot displacements (fault + fold + block vs. scaled fault : the aliasing problem)
    UX_temp = uzfault + 1.26.*uzfold_ds + 0.0.*uzfold_ts + 0.0.*uzblock + uzdetach;
    UZ_temp = (uzfault + uzdetach);
    scale = UZ_temp(:) \ UX_temp(:);
    UZ_temp = scale*UZ_temp;
    
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('vertical displacement aliasing');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-g', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-w', 'Linewidth', 3);
    plot(sx, UZ_temp, '--g', 'Linewidth', 3);    
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-g', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-w', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '--g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_z (all)', 'FontSize', 14);
    text(-22, -1.3-0.03, sprintf('u_z %2.1f(ramp+detachment)', scale), 'FontSize', 14);
    
    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
    %fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
    %set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
    plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
    plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    plot([0 -xdetach], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);

    
% Plot displacements (long term)
    UX_temp = uxlong;
    UZ_temp = uzlong;    
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('long term');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-w', 'Linewidth', 3);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);    
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-w', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, sprintf('u_z'), 'FontSize', 14);
    
    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
%     plot([0 -30], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);

    
    
    
% Plot interseismic velocities (long term)
    %% Messing around to look for smooth interseismic
    UX_temp = uxfault + 1.26.*uxfold_ds + 0.*uxfold_ts + 0.*uxblock + uxdetach;
    UZ_temp = uzfault + 1.26.*uzfold_ds + 0.*uzfold_ts + 0.*uzblock + uzdetach;
    %UX_temp = uxfault + uxfold_ds + uxfold_ts + uxblock + uxdetach;
    %UZ_temp = uzfault + uzfold_ds + uzfold_ts + uzblock + uzdetach;
    UX_total_all = UX_temp;
    UZ_total_all = UZ_temp;

    UX_temp = uxlong-1.0*UX_total_all;
    UZ_temp = uzlong-1.0*UZ_total_all;
    UX_new = UX_temp;
    UZ_new = UZ_temp;
 
   
    
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('long term - all');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-w', 'Linewidth', 3);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);    
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-w', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, sprintf('u_z'), 'FontSize', 14);
    
    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
%     plot([0 -30], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);
    
    
    
    
    
% Plot interseismic velocities (long term)
    UX_temp = uxlong-(uxfault+uxdetach);
    UZ_temp = uzlong-(uzfault+uzdetach);
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('long term - (ramp + detach)');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-w', 'Linewidth', 3);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);    
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-w', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, sprintf('u_z'), 'FontSize', 14);
    
    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
%     plot([0 -30], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);
    
    
    
% Complete disslocate vs classic solution
    uxclassic = uxblock - uxfault;
    uzclassic = uzblock - uzfault;
    
    % Unit block displacment
    idx                           = find(sx < faultx1);
    temp                          = R*[1 ; 0];
    uxblock                       = zeros(size(sx));
    uzblock                       = zeros(size(sx));
    uxblock(idx)                  = temp(1);
    uzblock(idx)                  = temp(2); 
%     uxblock(idx)                  = temp(1);
%     uzblock(idx)                  = 0; 
    uxblock = uxblock(:);
    uzblock = uzblock(:);
    % Classic block model
    uxclassic = (uxblock - uxfault)/uxblock(1); % have to renormalize
    uzclassic = (uzblock - uzfault)/uzblock(1);

    % Modified block model
%     uxclassic = (uxblock/uxblock(1) - uxfault); % have to renormalize
%     uzclassic = (uzblock/uzblock(1) - uzfault);

    
    
    % Plot interseismic velocities (long term)
    UX_temp = uxclassic;
    UZ_temp = uzclassic;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('long term - (ramp + detach)');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-w', 'Linewidth', 3);
    plot(sx, UZ_temp, '-g', 'Linewidth', 3);    
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-w', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, sprintf('u_z'), 'FontSize', 14);
    
    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
%     plot([0 -30], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);
    
    
% Classical and new horizontal interseismic
    % Plot interseismic velocities (long term)
    UX_temp = UX_new;
    UZ_temp = uxclassic;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('Interseismic, horizontal (new & classic)');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-b', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-w', 'Linewidth', 3);
    plot(sx, UZ_temp, '--b', 'Linewidth', 3);    
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-b', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-w', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '--b', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_x', 'FontSize', 14);
    text(-22, -1.3-0.03, sprintf('u_x'), 'FontSize', 14);
    
    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
%     plot([0 -30], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);


% Classical and new vertical interseismic
    UX_temp = UZ_new;
    UZ_temp = uzclassic;
    figure;
    ax1 = axes('Position', [0.2 0.4 0.6 0.4]); hold on;
    title('Interseismic, vertical (new & classic)');
    fh = fill([foldx2 faultx1 faultx1 foldx2], [-1.5 -1.5 1.5 1.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);
    plot(sx, UX_temp, '-k', 'Linewidth', 5);
    plot(sx, UX_temp, '-g', 'Linewidth', 3);
    plot(sx, UZ_temp, '-k', 'Linewidth', 5);
    plot(sx, UZ_temp, '-w', 'Linewidth', 3);
    plot(sx, UZ_temp, '--g', 'Linewidth', 3);    
    ylabel('u/u_0');
    set(gca, 'Xlim', [xmin xmax]);
    set(gca, 'XTick', []);
    set(gca, 'ylim', [-1.5 1.5])
    set(gca, 'ytick', [-1.5 -1.0 -0.5 0.0 0.5 1.0 1.5]);
    set(gca, 'ytickLabel', {'-1.5' '-1.0' '-0.5' '0.0' '0.5' '1.0' '1.5'});

    % Handmade legend
    plot([-23 -28], [-1 -1], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1 -1], '-g', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '-k', 'Linewidth', 5);
    plot([-23 -28], [-1.3 -1.3], '-w', 'Linewidth', 3);
    plot([-23 -28], [-1.3 -1.3], '--g', 'Linewidth', 3);
    text(-22, -1-0.03, 'u_z', 'FontSize', 14);
    text(-22, -1.3-0.03, sprintf('u_z'), 'FontSize', 14);
    
    % Cross section plot with active structures
    ax2 = axes('Position', [0.2 0.2 0.6 0.2]); hold on;

    % Fault and fold trace
    fh = fill([foldx2 faultx1 faultx1 foldx2], [12.5 12.5 -7.5 -7.5], 'y');
    set(fh, 'FaceColor', 0.75*[1 1 1], 'EdgeColor', 'none')

    % Fold limbs, if active
%     fh = fill([faultx1 faultx2 foldx2], [faultz1 faultz2 foldz2], 'r');
%     set(fh, 'EdgeColor', 'none');

    % Surface
    plot([xmin xmax], [0 0], '--k', 'Linewidth', 1, 'Color', 0.0*[1 1 1]);

    % Plot layer indicators to check layer thickness
    depvec = linspace(0, faultz2, nlay);
    latvec = linspace(foldx2, faultx1, nlay);
    for i = 2:numel(depvec)-1
       ytemp = depvec(i);
       xtemp = (ytemp+foldz2)/(-foldm) + foldx2; % Calculate intersection between horizontal and kink axis (hardwired #s that need to change)
       plot([xmin xtemp], [ytemp ytemp], '-k', 'Color', 0.5*[1 1 1]);
       plot([xtemp latvec(i)], [ytemp 0], '-k', 'Color', 0.5*[1 1 1]);
    end

    % Fault
    plot([faultx1 faultx2], [faultz1 faultz2], '-k', 'Linewidth', 5);
%     plot([faultx1 faultx2], [faultz1 faultz2], '-r', 'Linewidth', 3);

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);
%     plot([foldx1 foldx2], [foldz1 foldz2], '-r', 'Linewidth', 3);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
%     plot([0 -30], [faultz2 faultz2], '-r', 'Linewidth', 2);

    % Circles for cleaner looking intersections
    plot([faultx1 faultx2], [faultz1 faultz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot([foldx1 foldx2], [foldz1 foldz2], 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
    plot(-xdetach, faultz2, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

    set(gca, 'Xlim', [xmin xmax]);
    xlabel('x (km)');
    ylabel('z (km)');
    set(gca, 'Ydir', 'reverse');
    axis equal;
    set(gca, 'YTick', [0 5 10]);
    set(gca, 'ylim', [-7.5 12.5]);    
