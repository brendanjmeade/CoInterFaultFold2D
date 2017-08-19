close all; clear all;

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

% Calls to classical Okada for the horizontal case only.
xdetach                       = 8.6023; % Same length as ramp...at least for now
                                        %xdetach                       = 0; % Same length as ramp...at least for now

[uxdetach, ~, uzdetach]       = okada_plus_op(0, -5000, pi/2, 5, 0, 10000, xdetach, 0, -1, 0, sx, zeros(size(sx)), 0.25);

% Unit block displacment
idx                          = find(sx < faultx1);
temp                          = R*[1 ; 0];
uxblock                       = zeros(size(sx));
uzblock                       = zeros(size(sx));
uxblock(idx)                  = temp(1);
uzblock(idx)                  = temp(2);

% Total displacements
uxfault                       = uxfault(:);
uxblock                       = uxblock(:);
uzfault                       = uzfault(:);
uzblock                       = uzblock(:);

% Thick skinned case
%uxlong                        = uxblock;
%uxlong(sx < foldx2)           = 1;
%uzlong                        = uzblock;


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

    % Fold axis
    plot([foldx1 foldx2], [foldz1 foldz2], '-k', 'Linewidth', 5);

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);
    
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

    % Detachment
    plot([xmin 0], [faultz2 faultz2], '-k', 'Linewidth', 4);

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

    
