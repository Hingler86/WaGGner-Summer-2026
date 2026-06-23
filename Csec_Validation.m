%% Csec_Validation.m
% Acquisition time validation — crosshair (Csec) probe sweeps.
%
% Positions are assigned from Cmd_Delta_csec.txt sequentially
% (coordinate columns in data files are unreliable, ignored).
% ta=20s is used as the reference for error calculations.

clear; clc; close all;


%% USER SETTINGS

DATA_DIR = pwd;
CMD_PATH = "Cmd_Delta_csec.txt";

FILES = {
    "Csec_x500_P60_ta3_t1",   "ta = 3 s";
    "Csec_x500_P60_ta5_t1",   "ta = 5 s";
    "Csec_x500_P60_ta10_t1",  "ta = 10 s";
    "Csec_x500_P60_ta20_t1",  "ta = 20 s";
};

N_DATA_ROWS = 53;   % rows to use per file
REF_IDX     = 4;    % index of reference file (ta=20s)

X_POS_MM  = 500;
CMD_LEVEL = 60;

% Column indices (1-based)
COL_VMAG = 13;
COL_U = 10; COL_V = 11; COL_W = 12;

V_AMBIENT = 2.0;   % m/s — below this treated as out-of-jet


%% SHARED AXIS LIMITS (keeps all side-by-side plots comparable)

VMAG_YLIM  = [0, 10];    % m/s
ERR_YLIM   = [0, 30];    % % error (capped; dropouts excluded by dual in-jet mask)
POS_XLIM   = [-210, 210]; % mm


%% LOAD CMD PATH

cmd_raw = readmatrix(fullfile(DATA_DIR, CMD_PATH), ...
    "FileType","text","Delimiter","\t","NumHeaderLines",0);
cmd_Y = cmd_raw(:,1);
cmd_Z = cmd_raw(:,2);

y_rows = find(cmd_Z == 0 & cmd_Y ~= 0);
z_rows = find(cmd_Y == 0 & cmd_Z ~= 0);
Y_pos  = cmd_Y(y_rows);
Z_pos  = cmd_Z(z_rows);

fprintf("Cmd path: %d Y-sweep pts, %d Z-sweep pts\n", numel(Y_pos), numel(Z_pos));


%% LOAD ALL FILES

nFiles = size(FILES, 1);
D = struct([]);

for k = 1:nFiles
    raw = loadCsecFile(fullfile(DATA_DIR, FILES{k,1}), N_DATA_ROWS);
    D(k).label  = FILES{k,2};
    D(k).Vmag_Y = raw(y_rows, COL_VMAG);
    D(k).Vmag_Z = raw(z_rows, COL_VMAG);
    D(k).U_Y = raw(y_rows, COL_U); D(k).V_Y = raw(y_rows, COL_V); D(k).W_Y = raw(y_rows, COL_W);
    D(k).U_Z = raw(z_rows, COL_U); D(k).V_Z = raw(z_rows, COL_V); D(k).W_Z = raw(z_rows, COL_W);
    D(k).metrics = computeMetrics(Y_pos, D(k).Vmag_Y, Z_pos, D(k).Vmag_Z, V_AMBIENT);
    fprintf("Loaded %-30s  peak = %.3f m/s\n", FILES{k,1}, D(k).metrics.vmax);
end

tag    = sprintf("x=%d mm | %d%% cmd", X_POS_MM, CMD_LEVEL);
colors = lines(nFiles);
ref    = D(REF_IDX);


%% FIGURE 1 — Vmag profiles overlaid (shared axes)

figure("Name","Vmag Profiles — Acquisition Time Comparison", ...
    "Position",[50 50 1200 520]);

for col = 1:2
    subplot(1,2,col); hold on
    for k = 1:nFiles
        if col == 1
            plot(Y_pos, D(k).Vmag_Y, "-o", "Color",colors(k,:), ...
                "LineWidth",2,"MarkerSize",6,"MarkerFaceColor",colors(k,:), ...
                "DisplayName",D(k).label)
        else
            plot(Z_pos, D(k).Vmag_Z, "-o", "Color",colors(k,:), ...
                "LineWidth",2,"MarkerSize",6,"MarkerFaceColor",colors(k,:), ...
                "DisplayName",D(k).label)
        end
    end
    if col == 1
        xlabel("Y position (mm)"); title(sprintf("Y profile (Z=0) | %s", tag))
    else
        xlabel("Z position (mm)"); title(sprintf("Z profile (Y=0) | %s", tag))
    end
    ylabel("V_{mag} (m/s)")
    xlim(POS_XLIM); ylim(VMAG_YLIM)
    xline(0,"k:","LineWidth",0.8)
    legend("Location","northwest"); grid on; box on
end
sgtitle("Vmag profiles — acquisition time comparison","FontSize",13,"FontWeight","bold")


%% FIGURE 2 — Point-by-point % error vs reference (shared axes)

figure("Name","Point-by-point Error vs Reference", ...
    "Position",[50 620 1200 520]);

for col = 1:2
    subplot(1,2,col); hold on
    for k = 1:nFiles
        if k == REF_IDX; continue; end

        if col == 1
            test_v = D(k).Vmag_Y; ref_v = ref.Vmag_Y;
            pos = Y_pos;
        else
            test_v = D(k).Vmag_Z; ref_v = ref.Vmag_Z;
            pos = Z_pos;
        end

        % Only compute error where BOTH test and reference are in-jet
        in_jet = (ref_v > V_AMBIENT) & (test_v > V_AMBIENT);
        err = nan(size(pos));
        err(in_jet) = 100 * abs(test_v(in_jet) - ref_v(in_jet)) ./ ref_v(in_jet);

        plot(pos, err, "-o", "Color",colors(k,:), ...
            "LineWidth",2,"MarkerSize",6,"MarkerFaceColor",colors(k,:), ...
            "DisplayName",sprintf("%s vs %s", D(k).label, ref.label))
    end
    yline(5,  "k:","5%", "LineWidth",1.2,"LabelHorizontalAlignment","left")
    yline(10, "k--","10%","LineWidth",1.2,"LabelHorizontalAlignment","left")

    if col == 1
        xlabel("Y position (mm)")
        title(sprintf("Y profile error vs %s | %s", ref.label, tag))
    else
        xlabel("Z position (mm)")
        title(sprintf("Z profile error vs %s | %s", ref.label, tag))
    end
    ylabel("Error (%)")
    xlim(POS_XLIM); ylim(ERR_YLIM)
    xline(0,"k:","LineWidth",0.8)
    legend("Location","northeast"); grid on; box on
end
sgtitle(sprintf("Point-by-point error relative to reference: %s", ref.label), ...
    "FontSize",13,"FontWeight","bold")


%% FIGURE 3 — Velocity components per acquisition time (shared axes)

% Compute shared component axis limits
all_comp = [];
for k = 1:nFiles
    all_comp = [all_comp; D(k).U_Y; D(k).V_Y; D(k).W_Y; ...
                          D(k).U_Z; D(k).V_Z; D(k).W_Z]; %#ok<AGROW>
end
comp_lim = [min(all_comp)*1.1, max(all_comp)*1.1];

figure("Name","Velocity Components","Position",[50 50 1300 220*nFiles]);
comp_colors = [0.2 0.4 0.8; 0.8 0.2 0.2; 0.1 0.7 0.3];

for row = 1:nFiles
    d = D(row);
    subplot(nFiles,2,(row-1)*2+1); hold on
    plot(Y_pos, d.U_Y,"-","Color",comp_colors(1,:),"LineWidth",1.8,"DisplayName","U (axial)")
    plot(Y_pos, d.V_Y,"-","Color",comp_colors(2,:),"LineWidth",1.8,"DisplayName","V (lat. Y)")
    plot(Y_pos, d.W_Y,"-","Color",comp_colors(3,:),"LineWidth",1.8,"DisplayName","W (lat. Z)")
    xlabel("Y (mm)"); ylabel("Velocity (m/s)")
    title(sprintf("Y components | %s | %s", d.label, tag))
    xlim(POS_XLIM); ylim(comp_lim)
    legend("Location","northwest"); grid on; box on

    subplot(nFiles,2,(row-1)*2+2); hold on
    plot(Z_pos, d.U_Z,"-","Color",comp_colors(1,:),"LineWidth",1.8,"DisplayName","U (axial)")
    plot(Z_pos, d.V_Z,"-","Color",comp_colors(2,:),"LineWidth",1.8,"DisplayName","V (lat. Y)")
    plot(Z_pos, d.W_Z,"-","Color",comp_colors(3,:),"LineWidth",1.8,"DisplayName","W (lat. Z)")
    xlabel("Z (mm)"); ylabel("Velocity (m/s)")
    title(sprintf("Z components | %s | %s", d.label, tag))
    xlim(POS_XLIM); ylim(comp_lim)
    legend("Location","northwest"); grid on; box on
end
sgtitle("Velocity components per acquisition time","FontSize",13,"FontWeight","bold")


%% FIGURE 4 — Jet metrics bar chart (shared y-axes per metric)

figure("Name","Jet Metrics Summary","Position",[800 50 900 650]);

metric_keys  = {"vmax","fwhm_y","fwhm_z","yc","zc","std_y"};
metric_names = {"Peak V_{mag} (m/s)","FWHM Y (mm)","FWHM Z (mm)", ...
                "Centroid Y (mm)","Centroid Z (mm)","Spatial std Y (m/s)"};
labels = {D.label};

for m = 1:numel(metric_keys)
    vals = arrayfun(@(k) D(k).metrics.(metric_keys{m}), 1:nFiles);
    subplot(3,2,m)
    b = bar(vals,"FaceColor","flat");
    for k = 1:nFiles; b.CData(k,:) = colors(k,:); end
    set(gca,"XTickLabel",labels,"XTickLabelRotation",15)
    ylabel(metric_names{m}); title(metric_names{m})
    % Shared y-axis: pad by 20% above max
    ypad = max(abs(vals))*0.2;
    if ypad == 0; ypad = 0.5; end
    ylim([min(0, min(vals)-ypad), max(vals)+ypad])
    for b2 = 1:numel(vals)
        text(b2, vals(b2)+ypad*0.3, sprintf("%.2f",vals(b2)), ...
            "HorizontalAlignment","center","FontSize",9)
    end
    grid on; box on
end
sgtitle("Jet metrics vs acquisition time","FontSize",13,"FontWeight","bold")


%% CONSOLE SUMMARY

%fprintf("\n%s\n", repmat("=",1,65));
fprintf("ACQUISITION TIME VALIDATION — x=%dmm | %d%% cmd\n", X_POS_MM, CMD_LEVEL);
%fprintf("%s\n", repmat("=",1,65));
fprintf("%-12s %8s %8s %8s %8s %8s\n","Label","Vmax","Yc","Zc","FWHM_Y","FWHM_Z");
%fprintf("%s\n", repmat("-",1,65));
for k = 1:nFiles
    m = D(k).metrics;
    fprintf("%-12s %8.3f %8.2f %8.2f %8.1f %8.1f\n", ...
        D(k).label, m.vmax, m.yc, m.zc, m.fwhm_y, m.fwhm_z);
end

fprintf("\nIn-jet error vs reference (%s):\n", ref.label);
fprintf("%-12s %12s %12s %12s %12s\n","Label","Y mean","Y max","Z mean","Z max");
%fprintf("%s\n", repmat("-",1,65));
for k = 1:nFiles
    if k == REF_IDX; continue; end
    in_y = (ref.Vmag_Y > V_AMBIENT) & (D(k).Vmag_Y > V_AMBIENT);
    in_z = (ref.Vmag_Z > V_AMBIENT) & (D(k).Vmag_Z > V_AMBIENT);
    ey = 100*abs(D(k).Vmag_Y(in_y)-ref.Vmag_Y(in_y))./ref.Vmag_Y(in_y);
    ez = 100*abs(D(k).Vmag_Z(in_z)-ref.Vmag_Z(in_z))./ref.Vmag_Z(in_z);
    fprintf("%-12s %11.1f%% %11.1f%% %11.1f%% %11.1f%%\n", ...
        D(k).label, mean(ey), max(ey), mean(ez), max(ez));
end


%% LOCAL FUNCTIONS


function raw = loadCsecFile(filepath, nRows)
    lines = readlines(filepath);
    lines = lines(strtrim(lines) ~= "");
    lines = lines(1:min(nRows, numel(lines)));
    nCols = numel(split(lines(1), sprintf("\t")));
    data  = zeros(numel(lines), nCols);
    for k = 1:numel(lines)
        parts = split(lines(k), sprintf("\t"));
        for c = 1:min(numel(parts),nCols)
            s = strtrim(replace(string(parts(c)),",","."));
            v = str2double(s);
            if ~isnan(v); data(k,c) = v; end
        end
    end
    raw = data;
end

function m = computeMetrics(Y_pos, Vmag_Y, Z_pos, Vmag_Z, V_ambient)
    [m.vmax_y, iy] = max(Vmag_Y); [m.vmax_z, iz] = max(Vmag_Z);
    m.vmax = max(m.vmax_y, m.vmax_z);
    m.ymax = Y_pos(iy); m.zmax = Z_pos(iz);
    m.yc   = sum(Y_pos.*Vmag_Y)/sum(Vmag_Y);
    m.zc   = sum(Z_pos.*Vmag_Z)/sum(Vmag_Z);
    m.fwhm_y  = jetWidth(Y_pos, Vmag_Y, m.vmax_y/2);
    m.fwhm_z  = jetWidth(Z_pos, Vmag_Z, m.vmax_z/2);
    m.w1e_y   = jetWidth(Y_pos, Vmag_Y, m.vmax_y/exp(1));
    m.w1e_z   = jetWidth(Z_pos, Vmag_Z, m.vmax_z/exp(1));
    m.asym_y  = max(Vmag_Y(Y_pos>0)) / max(Vmag_Y(Y_pos<0));
    m.asym_z  = max(Vmag_Z(Z_pos>0)) / max(Vmag_Z(Z_pos<0));
    m.std_y   = std(Vmag_Y(Vmag_Y > V_ambient));
    m.std_z   = std(Vmag_Z(Vmag_Z > V_ambient));
end

function w = jetWidth(pos, vmag, threshold)
    f = vmag - threshold;
    edges = [];
    for k = 1:length(f)-1
        if f(k)*f(k+1) <= 0
            t = f(k)/(f(k)-f(k+1));
            edges(end+1) = pos(k)+t*(pos(k+1)-pos(k)); %#ok<AGROW>
        end
    end
    if numel(edges) >= 2; w = max(edges)-min(edges);
    else
        above = pos(vmag>=threshold);
        if numel(above)>=2; w = max(above)-min(above); else; w = NaN; end
    end
end