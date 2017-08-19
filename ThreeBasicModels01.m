close all; clear all;

% Plotting parameters
ms = 7; % marker size
lw = 2; % line width
fs = 14; % Font Size
EdgeColor = 'none';

% Declare station positions
n_stations                    = 1000;
xmin                          = -30;
xmax                          = 30;
sx                            = linspace(xmin, xmax, n_stations);
FontName                      = 'Helvectica';

% Primary fault surface
faultx1                       = 5;
faultz1                       = 0;
faultx2                       = 0;
faultz2                       = 5;
faultm                        = -(faultz2-faultz1)/(faultx2-faultx1);
gamma                         = abs(atand(-(faultz2-faultz1)/(faultx2-faultx1)));
R                             = [cosd(gamma) -sind(gamma) ; sind(gamma) cosd(gamma)]; % rotation matrix for 
bisectorAngle                 = (180-gamma)/2;
[uxfault, uzfault]            = DipTwoD(faultx1, faultz1, faultx2, faultz2, sx, 1, 0);

% Calls to classical Okada for the horizontal case only.
xdetach                       = 8.6023; % Same length as ramp...at least for now
[uxdetach, ~, uzdetach]       = okada_plus_op(0, -5000, pi/2, 5, deg2rad(135), 10000, 5*sqrt(2), 0, 0, 1, sx, zeros(size(sx)), 0.25);
kf1 = 1/cosd(45) - cosd(45);
kf2 = sind(45);

% Thick skinned case
idx                           = find(sx < faultx1);
temp                          = R*[sqrt(2) ; 0];
uxblock_thick                 = zeros(size(sx));
uzblock_thick                 = zeros(size(sx));
uxblock_thick(idx)            = temp(1);
uzblock_thick(idx)            = temp(2);

% convert all to column vectors
load('MakeMatU.mat');
bjm_ux_co = bjm_ux_co(:);
bjm_uz_co = bjm_uz_co(:);
bjm_ux_inter = bjm_ux_inter(:);
bjm_uz_inter = bjm_uz_inter(:);
sx = sx(:);
uxblock_thick = uxblock_thick(:);
uzblock_thick = uzblock_thick(:);
uxdetach = uxdetach(:);
uzdetach = uzdetach(:);
uxfault = uxfault(:);
uzfault = uzfault(:);


% Plot horizontal interseismic
figure('Position', [0 0 800 800]);
subplot(3, 2, 1); hold on;
plot([-30 30], [0 0], '-k', 'Color', 0.75*[1 1 1]);
plot([5 5], [-0.8 1.2], '-k', 'Color', 0.75*[1 1 1]);
plot(sx(:), uxblock_thick(:)-sqrt(2)*uxfault(:), '-g', 'LineWidth', lw);
plot(sx(:), uxblock_thick(:)-sqrt(2)*uxfault(:), '-b', 'LineWidth', lw);
plot(sx(:), uxblock_thick(:)-(kf1*uxfault(:) - kf2*uxdetach(:)), '-c', 'LineWidth', lw);
plot(sx(:), bjm_ux_inter(:), '-r', 'LineWidth', lw);
h1 = plot(sx(1:200:1000), uxblock_thick(1:200:1000)-sqrt(2)*uxfault(1:200:1000), 'gs', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h2 = plot(sx(50:200:1000), uxblock_thick(50:200:1000)-sqrt(2)*uxfault(50:200:1000), 'bs', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h3 = plot(sx(100:200:1000), uxblock_thick(100:200:1000)-(kf1*uxfault(100:200:1000) - kf2*uxdetach(100:200:1000)), 'cs', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h4 = plot(sx(150:200:1000), bjm_ux_inter(150:200:1000), 'rs', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
xlabel('x (km)')
ylabel('u/u_0')
set(gca, 'XLim', [-30 30]);
set(gca, 'YLim', [-0.8 1.2]);
set(gca, 'XTick', -30:15:30);
set(gca, 'YTick', -0.8:0.4:1.2);
ticks_format('%d', '%3.1f');
th = text(-28, -0.65, 'a) interseismic vertical', 'FontSize', fs, 'BackGroundColor', 'w');
lh = legend([h1 h2 h3 h4], 'model a', 'model b', 'model c', 'model d');
legend boxoff;
set(lh, 'Fontsize', 10);


% Plot vertical interseismic
subplot(3, 2, 2); hold on;
plot([-30 30], [0 0], '-k', 'Color', 0.75*[1 1 1]);
plot([5 5], [-0.8 1.2], '-k', 'Color', 0.75*[1 1 1]);
plot(sx(:), uzblock_thick(:)-sqrt(2)*uzfault(:), '-g', 'Linewidth', lw)
plot(sx(:), -sqrt(2)*uzfault(:), '-b', 'Linewidth', lw)
plot(sx(:), -(kf1*uzfault(:) - kf2*uzdetach(:)), '-c', 'Linewidth', lw)
plot(sx(:), bjm_uz_inter(:), '-r', 'LineWidth', lw);
h1 = plot(sx(1:200:1000), uzblock_thick(1:200:1000)-sqrt(2)*uzfault(1:200:1000), 'gs', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h2 = plot(sx(50:200:1000), -sqrt(2)*uzfault(50:200:1000), 'bs', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h3 = plot(sx(100:200:1000), -(kf1*uzfault(100:200:1000) - kf2*uzdetach(100:200:1000)), 'cs', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h4 = plot(sx(150:200:1000), bjm_uz_inter(150:200:1000), 'rs', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
xlabel('x (km)')
ylabel('u/u_0')
set(gca, 'XLim', [-30 30]);
set(gca, 'YLim', [-0.8 1.2]);
set(gca, 'XTick', -30:15:30);
set(gca, 'YTick', -0.8:0.4:1.2);
ticks_format('%d', '%3.1f');
th = text(-28, -0.65, 'b) interseismic vertical', 'FontSize', fs, 'BackGroundColor', 'w');
lh = legend([h1 h2 h3 h4], 'model a', 'model b', 'model c', 'model d');
legend boxoff;
set(lh, 'Fontsize', 10);


% Plot horizontal coseismic + interseismic
subplot(3, 2, 3); hold on;
plot([-30 30], [0 0], '-k', 'Color', 0.75*[1 1 1]);
plot([5 5], [-0.8 1.2], '-k', 'Color', 0.75*[1 1 1]);
plot(sx(:), sqrt(2)*uxfault(:), '-g', 'Linewidth', lw);
plot(sx(:), sqrt(2)*uxfault(:), '-b', 'Linewidth', lw);
plot(sx(:), sqrt(2)*uxfault(:), '-c', 'LineWidth', lw);
plot(sx(:), bjm_ux_co(:), '-r', 'LineWidth', lw);
h1 = plot(sx(1:200:1000), sqrt(2)*uxfault(1:200:1000), 'gs', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h2 = plot(sx(50:200:1000), sqrt(2)*uxfault(50:200:1000), 'bs', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h3 = plot(sx(100:200:1000), sqrt(2)*uxfault(100:200:1000), 'cs', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h4 = plot(sx(150:200:1000), bjm_ux_co(150:200:1000), 'rs', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
xlabel('x (km)')
ylabel('u/u_0')
set(gca, 'XLim', [-30 30]);
set(gca, 'YLim', [-0.8 1.2]);
set(gca, 'XTick', -30:15:30);
set(gca, 'YTick', -0.8:0.4:1.2);
ticks_format('%d', '%3.1f');
th = text(-28, -0.65, 'c) coseismic horizontal', 'FontSize', fs, 'BackGroundColor', 'w');
lh = legend([h1 h2 h3 h4], 'model a', 'model b', 'model c', 'model d');
legend boxoff;
set(lh, 'Fontsize', 10);


% Plot vertical coseismic + interseismic
subplot(3, 2, 4); hold on;
plot([-30 30], [0 0], '-k', 'Color', 0.75*[1 1 1]);
plot([5 5], [-0.8 1.2], '-k', 'Color', 0.75*[1 1 1]);
plot(sx(:), sqrt(2)*uzfault(:), '-g', 'Linewidth', lw)
plot(sx(:), sqrt(2)*uzfault(:), '-b')
plot(sx(:), sqrt(2)*uzfault(:), '-c', 'Linewidth', lw)
plot(sx(:), bjm_uz_co(:), '-r', 'LineWidth', lw);
h1 = plot(sx(1:200:1000), sqrt(2)*uzfault(1:200:1000), 'gs', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h2 = plot(sx(50:200:1000), sqrt(2)*uzfault(50:200:1000), 'bs', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h3 = plot(sx(100:200:1000), sqrt(2)*uzfault(100:200:1000), 'cs', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h4 = plot(sx(150:200:1000), bjm_uz_co(150:200:1000), 'rs', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
xlabel('x (km)');
ylabel('u/u_0');
set(gca, 'XLim', [-30 30]);
set(gca, 'YLim', [-0.8 1.2]);
set(gca, 'XTick', -30:15:30);
set(gca, 'YTick', -0.8:0.4:1.2);
ticks_format('%d', '%3.1f');
th = text(-28, -0.65, 'd) coseismic vertical', 'FontSize', fs, 'BackGroundColor', 'w');
lh = legend([h1 h2 h3 h4], 'model a', 'model b', 'model c', 'model d');
legend boxoff;
set(lh, 'Fontsize', 10);

% Plot horizontal coseismic + interseismic
subplot(3, 2, 5); hold on;
plot([-30 30], [0 0], '-k', 'Color', 0.75*[1 1 1]);
plot([5 5], [-0.8 1.2], '-k', 'Color', 0.75*[1 1 1]);
plot(sx(:), uxblock_thick(:)-sqrt(2)*uxfault(:) + sqrt(2)*uxfault(:), '-g', 'Linewidth', lw);
plot(sx(:), uxblock_thick(:)-sqrt(2)*uxfault(:) + sqrt(2)*uxfault(:), '-b', 'Linewidth', lw);
plot(sx(:), uxblock_thick(:)-(kf1*uxfault(:) - kf2*uxdetach(:)) + sqrt(2)*uxfault(:), '-c', 'LineWidth', lw);
plot(sx(:), bjm_ux_inter(:)+bjm_ux_co(:), '-r', 'LineWidth', lw);
h1 = plot(sx(1:200:1000), uxblock_thick(1:200:1000)-sqrt(2)*uxfault(1:200:1000) + sqrt(2)*uxfault(1:200:1000), 'gs', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h2 = plot(sx(50:200:1000), uxblock_thick(50:200:1000)-sqrt(2)*uxfault(50:200:1000) + sqrt(2)*uxfault(50:200:1000), 'bs', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h3 = plot(sx(100:200:1000), uxblock_thick(100:200:1000)-(kf1*uxfault(100:200:1000) - kf2*uxdetach(100:200:1000)) + sqrt(2)*uxfault(100:200:1000), 'cs', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h4 = plot(sx(150:200:1000), bjm_ux_inter(150:200:1000)+bjm_ux_co(150:200:1000), 'rs', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
xlabel('x (km)')
ylabel('u/u_0')
set(gca, 'XLim', [-30 30]);
set(gca, 'YLim', [-0.8 1.2]);
set(gca, 'XTick', -30:15:30);
set(gca, 'YTick', -0.8:0.4:1.2);
ticks_format('%d', '%3.1f');
th = text(-28, -0.65, 'e) (interseismic+coseismic) horizontal', 'FontSize', fs, 'BackGroundColor', 'w');
lh = legend([h1 h2 h3 h4], 'model a', 'model b', 'model c', 'model d');
legend boxoff;
set(lh, 'Fontsize', 10);

% Plot vertical coseismic + interseismic
subplot(3, 2, 6); hold on;
plot([-30 30], [0 0], '-k', 'Color', 0.75*[1 1 1]);
plot([5 5], [-0.8 1.2], '-k', 'Color', 0.75*[1 1 1]);
plot(sx(:), uzblock_thick(:)-sqrt(2)*uzfault(:) + sqrt(2)*uzfault(:), '-g', 'Linewidth', lw)
plot(sx(:), -sqrt(2)*uzfault(:) + sqrt(2)*uzfault(:), '-b')
plot(sx(:), -(kf1*uzfault(:) - kf2*uzdetach(:)) + sqrt(2)*uzfault(:), '-c', 'Linewidth', lw)
plot(sx(:), bjm_uz_inter(:)+bjm_uz_co(:), '-r', 'LineWidth', lw);
h1 = plot(sx(1:200:1000), uzblock_thick(1:200:1000)-sqrt(2)*uzfault(1:200:1000) + sqrt(2)*uzfault(1:200:1000), 'gs', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h2 = plot(sx(50:200:1000), uzblock_thick(50:200:1000)-sqrt(2)*uzfault(50:200:1000) + sqrt(2)*uzfault(50:200:1000), 'bs', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h3 = plot(sx(100:200:1000), -(kf1*uzfault(100:200:1000) - kf2*uzdetach(100:200:1000)) + sqrt(2)*uzfault(100:200:1000), 'cs', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
h4 = plot(sx(150:200:1000), bjm_uz_inter(150:200:1000)+bjm_uz_co(150:200:1000), 'rs', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', EdgeColor, 'MarkerSize', ms);
xlabel('x (km)');
ylabel('u/u_0');
set(gca, 'XLim', [-30 30]);
set(gca, 'YLim', [-0.8 1.2]);
set(gca, 'XTick', -30:15:30);
set(gca, 'YTick', -0.8:0.4:1.2);
ticks_format('%d', '%3.1f');
th = text(-28, -0.65, 'f) (interseismic+coseismic) vertical', 'FontSize', fs, 'BackGroundColor', 'w');
lh = legend([h1 h2 h3 h4], 'model a', 'model b', 'model c', 'model d');
legend boxoff;
set(lh, 'Fontsize', 10);



