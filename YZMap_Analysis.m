%% YZMap_Analysis.m
% Jet characterization from YZ mapping (five-hole probe, auto mode)
% Supports both cross-pattern and full snake-grid scan paths.
%
% Input file format (tab-separated, no header):
%   Col 1:  Pos Y (mm)
%   Col 2:  Pos Z (mm)
%   Col 3-7: P1-P5 - Pamb (hPa)
%   Col 8:  Pamb (hPa)
%   Col 9:  Tamb (°C)
%   Col 10: U  (m/s)  — axial
%   Col 11: V  (m/s)  — lateral Y
%   Col 12: W  (m/s)  — lateral Z
%   Col 13: Norme V / Vmag (m/s)
%   Col 14: Alpha (°)
%   Col 15-19: Beta, Theta, Phi, Psi (°)
%   Col 20: Mach
%   Col 21: Ptot abs (hPa)
%   Col 22: Ps abs (hPa)
%   Col 23: Angle Sonde (deg)

clear; clc; close all;


%% USER SETTINGS

DATA_DIR     = pwd;
FILE_PATTERN = "YZMAP_x100_P60_t1";   % filename or wildcard

X_POS_MM     = 100;    % downstream distance (mm) — for titles/labels
CMD_LEVEL    = 30;     % command level (%) — for titles/labels

% Column indices (1-based)
COL_Y    = 1;
COL_Z    = 2;
COL_U    = 10;   % axial velocity
COL_V    = 11;   % lateral Y velocity
COL_W    = 12;   % lateral Z velocity
COL_VMAG = 13;   % velocity magnitude

% Ambient threshold (m/s) — points below this are considered out-of-jet
V_AMBIENT = 2.5;

% Interpolation grid resolution for contour plots (mm)
GRID_RES = 3;


%% LOAD FILE

files = dir(fullfile(DATA_DIR, FILE_PATTERN));
if isempty(files)
    error("No files found matching: %s", FILE_PATTERN);
end

filepath = fullfile(files(1).folder, files(1).name);
raw = table2array(readYZMapFile(filepath));

fprintf("Loaded: %s  (%d measurement points)\n", files(1).name, size(raw,1));

Y    = raw(:, COL_Y);
Z    = raw(:, COL_Z);
U    = raw(:, COL_U);
V    = raw(:, COL_V);
W    = raw(:, COL_W);
Vmag = raw(:, COL_VMAG);


%% EXTRACT CARDINAL PROFILES (Z=0 row and Y=0 column)
%% These always exist in both cross and snake grids.

tol = 1e-6;  % floating point tolerance for position matching

% Y profile: all points where Z == 0, sorted by Y
y_mask = abs(Z) < tol;
[Y_prof, si] = sort(Y(y_mask));
Vmag_Yp = Vmag(y_mask); Vmag_Yp = Vmag_Yp(si);
U_Yp    = U(y_mask);    U_Yp    = U_Yp(si);
V_Yp    = V(y_mask);    V_Yp    = V_Yp(si);
W_Yp    = W(y_mask);    W_Yp    = W_Yp(si);

% Z profile: all points where Y == 0, sorted by Z
z_mask = abs(Y) < tol;
[Z_prof, si] = sort(Z(z_mask));
Vmag_Zp = Vmag(z_mask); Vmag_Zp = Vmag_Zp(si);
U_Zp    = U(z_mask);    U_Zp    = U_Zp(si);
V_Zp    = V(z_mask);    V_Zp    = V_Zp(si);
W_Zp    = W(z_mask);    W_Zp    = W_Zp(si);

fprintf("Y profile (Z=0): %d points\n", numel(Y_prof));
fprintf("Z profile (Y=0): %d points\n", numel(Z_prof));


%% BUILD INTERPOLATED 2D FIELD

Y_range = min(Y):GRID_RES:max(Y);
Z_range = min(Z):GRID_RES:max(Z);
[YG, ZG] = meshgrid(Y_range, Z_range);

F_vmag = scatteredInterpolant(Y, Z, Vmag, "natural", "nearest");
F_U    = scatteredInterpolant(Y, Z, U,    "natural", "nearest");
F_V    = scatteredInterpolant(Y, Z, V,    "natural", "nearest");
F_W    = scatteredInterpolant(Y, Z, W,    "natural", "nearest");

VG = F_vmag(YG, ZG);
UG = F_U(YG, ZG);
VvG = F_V(YG, ZG);
WG = F_W(YG, ZG);


%% JET METRICS

metrics = computeJetMetrics(Y_prof, Vmag_Yp, Z_prof, Vmag_Zp, Y, Z, Vmag, V_AMBIENT);

tag = sprintf("x=%d mm | %d%% cmd", X_POS_MM, CMD_LEVEL);

fprintf("\n========== JET CHARACTERIZATION ==========\n");
fprintf("File       : %s\n",    files(1).name);
fprintf("X position : %d mm\n", X_POS_MM);
fprintf("Command    : %d%%\n",  CMD_LEVEL);
fprintf("-------------------------------------------\n");
fprintf("Overall max Vmag   : %.3f m/s at (Y=%.0f, Z=%.0f) mm\n", ...
    metrics.Vmax_2D, metrics.Ymax_2D, metrics.Zmax_2D);
fprintf("Max along Y profile: %.3f m/s at Y = %.0f mm\n", metrics.Vmax_Y, metrics.Ymax);
fprintf("Max along Z profile: %.3f m/s at Z = %.0f mm\n", metrics.Vmax_Z, metrics.Zmax);
fprintf("Jet centre (velocity-weighted centroid, full 2D):\n");
fprintf("   Y_c = %.2f mm,  Z_c = %.2f mm\n", metrics.Y_center, metrics.Z_center);
fprintf("FWHM jet width (from cardinal profiles):\n");
fprintf("   Y direction: %.1f mm\n", metrics.FWHM_Y);
fprintf("   Z direction: %.1f mm\n", metrics.FWHM_Z);
fprintf("1/e jet width (from cardinal profiles):\n");
fprintf("   Y direction: %.1f mm\n", metrics.width1e_Y);
fprintf("   Z direction: %.1f mm\n", metrics.width1e_Z);
fprintf("Asymmetry (Vmax pos. / Vmax neg. half):\n");
fprintf("   Y axis: %.3f\n", metrics.asym_Y);
fprintf("   Z axis: %.3f\n", metrics.asym_Z);
fprintf("Spatial std dev (in-jet points, full 2D): %.3f m/s\n", metrics.std_2D);
fprintf("===========================================\n\n");


%% FIGURE 1 — Vmag cardinal profiles

figure("Name","Vmag Profiles","Position",[50 50 1100 480]);

subplot(1,2,1)
plot(Y_prof, Vmag_Yp, "b-o","LineWidth",1.8,"MarkerSize",6,"MarkerFaceColor","b")
hold on
xline(metrics.Y_center, "r--", sprintf("Y_c = %.1f mm", metrics.Y_center), ...
    "LineWidth",1.2,"LabelVerticalAlignment","bottom")
yline(metrics.Vmax_Y/2,      ":","FWHM level",  "Color",[0.5 0.5 0.5],"LineWidth",1)
yline(metrics.Vmax_Y/exp(1), "--","1/e level",   "Color",[0.7 0.4 0.0],"LineWidth",1)
xlabel("Y position (mm)"); ylabel("V_{mag} (m/s)")
title(sprintf("Y profile (Z=0) | %s", tag))
legend("V_{mag}","Jet centre","FWHM level","1/e level","Location","northwest")
grid on; box on

subplot(1,2,2)
plot(Z_prof, Vmag_Zp, "r-o","LineWidth",1.8,"MarkerSize",6,"MarkerFaceColor","r")
hold on
xline(metrics.Z_center, "b--", sprintf("Z_c = %.1f mm", metrics.Z_center), ...
    "LineWidth",1.2,"LabelVerticalAlignment","bottom")
yline(metrics.Vmax_Z/2,      ":","FWHM level",  "Color",[0.5 0.5 0.5],"LineWidth",1)
yline(metrics.Vmax_Z/exp(1), "--","1/e level",   "Color",[0.7 0.4 0.0],"LineWidth",1)
xlabel("Z position (mm)"); ylabel("V_{mag} (m/s)")
title(sprintf("Z profile (Y=0) | %s", tag))
legend("V_{mag}","Jet centre","FWHM level","1/e level","Location","northwest")
grid on; box on


%% FIGURE 2 — Velocity components along cardinal profiles

figure("Name","Velocity Components","Position",[50 580 1100 480]);

subplot(1,2,1)
plot(Y_prof, U_Yp, "b-o","LineWidth",1.5,"MarkerSize",5,"DisplayName","U (axial)")
hold on
plot(Y_prof, V_Yp, "r-s","LineWidth",1.5,"MarkerSize",5,"DisplayName","V (lateral Y)")
plot(Y_prof, W_Yp, "g-^","LineWidth",1.5,"MarkerSize",5,"DisplayName","W (lateral Z)")
xlabel("Y position (mm)"); ylabel("Velocity (m/s)")
title(sprintf("Velocity components — Y profile | %s", tag))
legend; grid on; box on

subplot(1,2,2)
plot(Z_prof, U_Zp, "b-o","LineWidth",1.5,"MarkerSize",5,"DisplayName","U (axial)")
hold on
plot(Z_prof, V_Zp, "r-s","LineWidth",1.5,"MarkerSize",5,"DisplayName","V (lateral Y)")
plot(Z_prof, W_Zp, "g-^","LineWidth",1.5,"MarkerSize",5,"DisplayName","W (lateral Z)")
xlabel("Z position (mm)"); ylabel("Velocity (m/s)")
title(sprintf("Velocity components — Z profile | %s", tag))
legend; grid on; box on


%% FIGURE 3 — 2D Vmag contour map

figure("Name","2D Vmag Map","Position",[200 100 750 640]);

contourf(YG, ZG, VG, 25, "LineColor","none")
hold on

% FWHM contour — use overall 2D peak
half_max = metrics.Vmax_2D / 2;
contour(YG, ZG, VG, [half_max half_max], "k--","LineWidth",1.8,"DisplayName","FWHM contour")

% Jet centroid
plot(metrics.Y_center, metrics.Z_center, "k+","MarkerSize",16, ...
    "LineWidth",2.5,"DisplayName","Jet centre")

% Measurement points
plot(Y, Z, "w.","MarkerSize",5,"DisplayName","Measured points")

colorbar
colormap("jet")
clim([0, max(Vmag)*1.05])
xlabel("Y (mm)"); ylabel("Z (mm)")
title(sprintf("V_{mag} (m/s) — 2D map | %s", tag))
legend("Location","northeast")
axis equal; grid on; box on


%% FIGURE 4 — 2D U, V, W component maps (3-panel)

figure("Name","2D Velocity Components","Position",[200 780 1200 440]);

comp_data  = {UG,   VvG,  WG};
comp_names = {"U (axial, m/s)", "V (lateral Y, m/s)", "W (lateral Z, m/s)"};

for k = 1:3
    subplot(1,3,k)
    contourf(YG, ZG, comp_data{k}, 20, "LineColor","none")
    hold on
    plot(Y, Z, "w.","MarkerSize",4)
    plot(metrics.Y_center, metrics.Z_center, "k+","MarkerSize",12,"LineWidth",2)
    colorbar; colormap("jet")
    xlabel("Y (mm)"); ylabel("Z (mm)")
    title(sprintf("%s | %s", comp_names{k}, tag))
    axis equal; grid on; box on
end


%% FIGURE 5 — Jet width summary

figure("Name","Jet Width Summary","Position",[1000 100 560 420]);

categories = ["FWHM Y", "FWHM Z", "1/e Y", "1/e Z"];
widths_mm  = [metrics.FWHM_Y, metrics.FWHM_Z, metrics.width1e_Y, metrics.width1e_Z];

b = bar(widths_mm, "FaceColor","flat");
b.CData = [0.2 0.4 0.8; 0.8 0.2 0.2; 0.4 0.6 1.0; 1.0 0.5 0.5];
set(gca,"XTickLabel", categories)
ylabel("Jet width (mm)")
title(sprintf("Jet width summary | %s", tag))
grid on; box on

for k = 1:numel(widths_mm)
    if ~isnan(widths_mm(k))
        text(k, widths_mm(k)+0.5, sprintf("%.0f mm", widths_mm(k)), ...
            "HorizontalAlignment","center","FontSize",10)
    end
end


%% LOCAL FUNCTIONS


function raw = readYZMapFile(filepath)
%READYZMAPFILE  Read auto-mode five-hole probe output (no header, comma decimals).
    lines = readlines(filepath);
    lines = lines(strtrim(lines) ~= "");
    nCols = numel(split(lines(1), sprintf("\t")));
    data  = zeros(numel(lines), nCols);
    for k = 1:numel(lines)
        parts = split(lines(k), sprintf("\t"));
        for c = 1:numel(parts)
            s = strtrim(replace(string(parts(c)), ",", "."));
            v = str2double(s);
            if isnan(v); v = 0; end
            data(k,c) = v;
        end
    end
    raw = array2table(data);
end


function metrics = computeJetMetrics(Y_prof, Vmag_Yp, Z_prof, Vmag_Zp, ...
                                      Y_all, Z_all, Vmag_all, V_ambient)
%COMPUTEJETMETRICS  Extract jet characterization quantities from profiles + full 2D data.

    %% Profile peaks
    [metrics.Vmax_Y, iY] = max(Vmag_Yp);
    [metrics.Vmax_Z, iZ] = max(Vmag_Zp);
    metrics.Ymax = Y_prof(iY);
    metrics.Zmax = Z_prof(iZ);

    %% Overall 2D peak
    [metrics.Vmax_2D, i2D] = max(Vmag_all);
    metrics.Ymax_2D = Y_all(i2D);
    metrics.Zmax_2D = Z_all(i2D);

    %% 2D velocity-weighted centroid
    metrics.Y_center = sum(Y_all .* Vmag_all) / sum(Vmag_all);
    metrics.Z_center = sum(Z_all .* Vmag_all) / sum(Vmag_all);

    %% FWHM (from cardinal profiles)
    metrics.FWHM_Y    = jetWidth(Y_prof, Vmag_Yp, metrics.Vmax_Y/2);
    metrics.FWHM_Z    = jetWidth(Z_prof, Vmag_Zp, metrics.Vmax_Z/2);

    %% 1/e width (from cardinal profiles)
    metrics.width1e_Y = jetWidth(Y_prof, Vmag_Yp, metrics.Vmax_Y/exp(1));
    metrics.width1e_Z = jetWidth(Z_prof, Vmag_Zp, metrics.Vmax_Z/exp(1));

    %% Asymmetry: peak velocity on positive vs negative half of each profile
    metrics.asym_Y = max(Vmag_Yp(Y_prof > 0)) / max(Vmag_Yp(Y_prof < 0));
    metrics.asym_Z = max(Vmag_Zp(Z_prof > 0)) / max(Vmag_Zp(Z_prof < 0));

    %% Spatial std dev across all in-jet 2D points
    in_jet = Vmag_all > V_ambient;
    metrics.std_2D = std(Vmag_all(in_jet));
end


function w = jetWidth(pos, vmag, threshold)
%JETWIDTH  Width of region where vmag >= threshold, via linear interpolation.
    above = vmag >= threshold;
    if sum(above) < 2
        w = NaN;
        return
    end
    f = vmag - threshold;
    edges = [];
    for k = 1:length(f)-1
        if f(k)*f(k+1) <= 0
            t = f(k) / (f(k) - f(k+1));
            edges(end+1) = pos(k) + t*(pos(k+1) - pos(k)); %#ok<AGROW>
        end
    end
    if numel(edges) < 2
        w = max(pos(above)) - min(pos(above));
    else
        w = max(edges) - min(edges);
    end
end