%% YZMap_Comparison.m
% Comparison figures across all X positions and command levels
% for presentation to research professor.

clear; clc; close all;


%% USER SETTINGS

DATA_DIR = pwd;

files = struct( ...
    'path',  { ...
        "YZMAP_x100_P30_t1", ...
        "YZMAP_x100_P90_t2", ...
        "YZMAP_x400_P30_t1", ...
        "YZMAP_x400_P90_t3", ...
        "YZMAP_x600_P30_t1", ...
        "YZMAP_x600_P90_t2"  ...
    }, ...
    'x',   {100, 100, 400, 400, 600, 600}, ...
    'pct', { 30,  90,  30,  90,  30,  90}  ...
);

COL_Y=1; COL_Z=2; COL_U=10; COL_V=11; COL_W=12; COL_VMAG=13;
V_AMBIENT = 2.0;
GRID_RES  = 3;  % mm interpolation resolution

% Plot style
colors30 = [0.20 0.45 0.80];   % blue family for 30%
colors90 = [0.85 0.20 0.20];   % red family for 90%
xMarkers = {'o','s','^'};       % marker per X position
xVals    = [100, 400, 600];


%% LOAD ALL DATA

D = struct([]);
for k = 1:numel(files)
    raw  = table2array(readYZMapFile(fullfile(DATA_DIR, files(k).path)));
    D(k).x    = files(k).x;
    D(k).pct  = files(k).pct;
    D(k).Y    = raw(:,COL_Y);
    D(k).Z    = raw(:,COL_Z);
    D(k).U    = raw(:,COL_U);
    D(k).V    = raw(:,COL_V);
    D(k).W    = raw(:,COL_W);
    D(k).Vmag = raw(:,COL_VMAG);
    D(k).metrics = computeMetrics(D(k), V_AMBIENT);
    fprintf("Loaded x=%d P=%d%%  peak=%.2f m/s\n", D(k).x, D(k).pct, D(k).metrics.vmax);
end


%% FIGURE 1 — 2D Vmag contour maps, 2x3 grid

fig1 = figure("Name","2D Vmag Maps — All Conditions","Position",[50 50 1400 800]);

plotOrder = [1 3 5; 2 4 6];  % row1=30%, row2=90%, cols=x100,x400,x600
pctLabels = ["30%","90%"];

for row = 1:2
    pct = [30 90];
    for col = 1:3
        idx = plotOrder(row,col);
        d   = D(idx);
        m   = d.metrics;

        % Interpolate
        Yr = min(d.Y):GRID_RES:max(d.Y);
        Zr = min(d.Z):GRID_RES:max(d.Z);
        [YG,ZG] = meshgrid(Yr,Zr);
        F  = scatteredInterpolant(d.Y, d.Z, d.Vmag, "natural","nearest");
        VG = F(YG,ZG);

        subplot(2,3,(row-1)*3+col)
        contourf(YG,ZG,VG,20,"LineColor","none")
        hold on
        contour(YG,ZG,VG,[m.vmax/2 m.vmax/2],"k--","LineWidth",1.5)
        plot(m.yc, m.zc, "k+","MarkerSize",14,"LineWidth",2.5)
        colorbar
        colormap("jet")
        clim([0 max([D(row*3-2).metrics.vmax, D(row*3-1).metrics.vmax, D(row*3).metrics.vmax])*1.05])
        xlabel("Y (mm)"); ylabel("Z (mm)")
        title(sprintf("x=%d mm | %d%% | peak=%.1f m/s", d.x, d.pct, m.vmax))
        axis equal; grid on; box on
    end
end
sgtitle("V_{mag} (m/s) — 2D maps | Top: 30% command | Bottom: 90% command", ...
    "FontSize",13,"FontWeight","bold")


%% FIGURE 2 — Y and Z profiles overlaid across X positions

fig2 = figure("Name","Cardinal Profiles — X Position Comparison","Position",[50 900 1400 700]);

for row = 1:2
    pct = [30 90]; p = pct(row);
    col_base = if_val(p==30, colors30, colors90);

    % Y profiles
    subplot(2,2,(row-1)*2+1)
    hold on
    for ki = 1:3
        idx = find([D.x]==xVals(ki) & [D.pct]==p);
        d = D(idx); m = d.metrics;
        lw = 2.5 - ki*0.3;
        alpha_shade = 1.0 - (ki-1)*0.25;
        c = col_base * alpha_shade + (1-alpha_shade)*[1 1 1];
        plot(d.metrics.Yp, d.metrics.Vp_y, ...
            "Color",c,"LineWidth",lw,"Marker",xMarkers{ki}, ...
            "MarkerSize",6,"MarkerFaceColor",c, ...
            "DisplayName",sprintf("x=%d mm (%.1f m/s)",xVals(ki),m.vmax))
    end
    xlabel("Y position (mm)"); ylabel("V_{mag} (m/s)")
    title(sprintf("Y profile (Z=0) | %d%% command", p))
    legend("Location","northwest"); grid on; box on
    xline(0,"k:","LineWidth",0.8)

    % Z profiles
    subplot(2,2,(row-1)*2+2)
    hold on
    for ki = 1:3
        idx = find([D.x]==xVals(ki) & [D.pct]==p);
        d = D(idx); m = d.metrics;
        lw = 2.5 - ki*0.3;
        alpha_shade = 1.0 - (ki-1)*0.25;
        c = col_base * alpha_shade + (1-alpha_shade)*[1 1 1];
        plot(d.metrics.Zp, d.metrics.Vp_z, ...
            "Color",c,"LineWidth",lw,"Marker",xMarkers{ki}, ...
            "MarkerSize",6,"MarkerFaceColor",c, ...
            "DisplayName",sprintf("x=%d mm (%.1f m/s)",xVals(ki),m.vmax))
    end
    xlabel("Z position (mm)"); ylabel("V_{mag} (m/s)")
    title(sprintf("Z profile (Y=0) | %d%% command", p))
    legend("Location","northwest"); grid on; box on
    xline(0,"k:","LineWidth",0.8)
end
sgtitle("Cardinal velocity profiles — downstream evolution", ...
    "FontSize",13,"FontWeight","bold")


%% FIGURE 3 — Peak velocity decay with downstream distance

fig3 = figure("Name","Peak Velocity Decay","Position",[500 50 750 550]);

hold on
vmax30 = arrayfun(@(k) D(k).metrics.vmax, find([D.pct]==30));
vmax90 = arrayfun(@(k) D(k).metrics.vmax, find([D.pct]==90));

plot(xVals, vmax30, "b-o","LineWidth",2.5,"MarkerSize",10, ...
    "MarkerFaceColor","b","DisplayName","30% command")
plot(xVals, vmax90, "r-s","LineWidth",2.5,"MarkerSize",10, ...
    "MarkerFaceColor","r","DisplayName","90% command")

% Annotate decay %
for ki = 2:3
    pct_decay30 = (1 - vmax30(ki)/vmax30(1))*100;
    pct_decay90 = (1 - vmax90(ki)/vmax90(1))*100;
    text(xVals(ki), vmax30(ki)-0.5, sprintf("−%.0f%%",pct_decay30), ...
        "HorizontalAlignment","center","Color","b","FontSize",10)
    text(xVals(ki), vmax90(ki)+0.8, sprintf("−%.0f%%",pct_decay90), ...
        "HorizontalAlignment","center","Color","r","FontSize",10)
end

xlabel("Downstream distance X (mm)")
ylabel("Peak V_{mag} (m/s)")
title("Peak velocity decay with downstream distance")
legend("Location","northeast"); grid on; box on
xlim([50 650]); xticks(xVals)


%% FIGURE 4 — 30% vs 90% profile comparison at each X position

fig4 = figure("Name","30% vs 90% Command Comparison","Position",[500 650 1200 650]);

for col = 1:3
    x = xVals(col);
    d30 = D(find([D.x]==x & [D.pct]==30));
    d90 = D(find([D.x]==x & [D.pct]==90));

    subplot(1,3,col)
    hold on

    % Normalize both to their own peak for shape comparison
    yyaxis left
    plot(d30.metrics.Yp, d30.metrics.Vp_y, "b-o","LineWidth",2,"MarkerSize",5, ...
        "MarkerFaceColor","b","DisplayName","30%")
    plot(d90.metrics.Yp, d90.metrics.Vp_y, "r-s","LineWidth",2,"MarkerSize",5, ...
        "MarkerFaceColor","r","DisplayName","90%")
    ylabel("V_{mag} (m/s)")
    ylim([0, max(d90.metrics.Vp_y)*1.15])

    xlabel("Y position (mm)")
    title(sprintf("Y profile (Z=0) | x=%d mm", x))
    legend("Location","northwest"); grid on; box on
    xline(0,"k:","LineWidth",0.8)
end
sgtitle("30% vs 90% command — Y profile shape comparison at each X position", ...
    "FontSize",13,"FontWeight","bold")


%% FIGURE 5 — Jet centroid drift

fig5 = figure("Name","Jet Centroid Position","Position",[1280 50 600 500]);

subplot(2,1,1)
hold on
yc30 = arrayfun(@(k) D(k).metrics.yc, find([D.pct]==30));
yc90 = arrayfun(@(k) D(k).metrics.yc, find([D.pct]==90));
plot(xVals, yc30, "b-o","LineWidth",2,"MarkerSize",8,"MarkerFaceColor","b","DisplayName","30%")
plot(xVals, yc90, "r-s","LineWidth",2,"MarkerSize",8,"MarkerFaceColor","r","DisplayName","90%")
yline(0,"k--","LineWidth",1,"Label","Geometric centre")
xlabel("X (mm)"); ylabel("Y centroid (mm)")
title("Y centroid position vs downstream distance")
legend; grid on; box on; xlim([50 650]); xticks(xVals)

subplot(2,1,2)
hold on
zc30 = arrayfun(@(k) D(k).metrics.zc, find([D.pct]==30));
zc90 = arrayfun(@(k) D(k).metrics.zc, find([D.pct]==90));
plot(xVals, zc30, "b-o","LineWidth",2,"MarkerSize",8,"MarkerFaceColor","b","DisplayName","30%")
plot(xVals, zc90, "r-s","LineWidth",2,"MarkerSize",8,"MarkerFaceColor","r","DisplayName","90%")
yline(0,"k--","LineWidth",1,"Label","Geometric centre")
xlabel("X (mm)"); ylabel("Z centroid (mm)")
title("Z centroid position vs downstream distance")
legend; grid on; box on; xlim([50 650]); xticks(xVals)
sgtitle("Jet centroid drift", "FontSize",13,"FontWeight","bold")


%% CONSOLE SUMMARY

fprintf("\n%s\n", repmat("=",1,70));
fprintf("PRESENTATION SUMMARY\n");
fprintf("%s\n", repmat("=",1,70));
for p = [30 90]
    fprintf("\n--- %d%% command ---\n",p);
    fprintf("  %-20s %8s %8s %8s\n","Metric","x=100","x=400","x=600");
    fprintf("  %s\n",repmat("-",1,48));
    idxs = find([D.pct]==p);
    metrics_list = { ...
        "Peak Vmag (m/s)",  "vmax"; ...
        "Y centroid (mm)",  "yc";   ...
        "Z centroid (mm)",  "zc";   ...
        "Asymmetry Y",      "asym_y"; ...
        "Asymmetry Z",      "asym_z"; ...
        "Spatial std (m/s)","std"};
    for m = 1:size(metrics_list,1)
        label = metrics_list{m,1};
        key   = metrics_list{m,2};
        vals  = arrayfun(@(k) D(k).metrics.(key), idxs);
        fprintf("  %-20s %8.3f %8.3f %8.3f\n", label, vals(1), vals(2), vals(3));
    end
end


%% LOCAL FUNCTIONS

function m = computeMetrics(d, V_ambient)
    tol = 1e-6;
    ym = abs(d.Z) < tol; zm = abs(d.Y) < tol;
    [Yp,si] = sort(d.Y(ym)); Vp_y = d.Vmag(ym); Vp_y = Vp_y(si);
    [Zp,si] = sort(d.Z(zm)); Vp_z = d.Vmag(zm); Vp_z = Vp_z(si);
    [vmax,i2d] = max(d.Vmag);
    m.vmax   = vmax;
    m.ymax   = d.Y(i2d); m.zmax = d.Z(i2d);
    m.yc     = sum(d.Y .* d.Vmag) / sum(d.Vmag);
    m.zc     = sum(d.Z .* d.Vmag) / sum(d.Vmag);
    m.vmax_y = max(Vp_y); m.vmax_z = max(Vp_z);
    m.asym_y = max(Vp_y(Yp>0)) / max(Vp_y(Yp<0));
    m.asym_z = max(Vp_z(Zp>0)) / max(Vp_z(Zp<0));
    m.std    = std(d.Vmag(d.Vmag > V_ambient));
    m.Yp = Yp; m.Vp_y = Vp_y;
    m.Zp = Zp; m.Vp_z = Vp_z;
end

function v = if_val(cond, a, b)
    if cond; v = a; else; v = b; end
end

function raw = readYZMapFile(filepath)
    lines = readlines(filepath);
    lines = lines(strtrim(lines) ~= "");
    nCols = numel(split(lines(1), sprintf("\t")));
    data  = zeros(numel(lines), nCols);
    for k = 1:numel(lines)
        parts = split(lines(k), sprintf("\t"));
        for c = 1:min(numel(parts),nCols)
            s = strtrim(replace(string(parts(c)), ",", "."));
            v = str2double(s);
            if ~isnan(v); data(k,c) = v; end
        end
    end
    raw = array2table(data);
end