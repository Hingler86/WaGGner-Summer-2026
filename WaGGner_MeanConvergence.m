%% WaGGner_MeanC.m
% Five-hole probe mean convergence analysis
% Simple version: individual figures + grouped power figures

clear; clc; close all;

%% ================= USER SETTINGS =================

DATA_DIR = pwd;
FILE_PATTERN = "MeanConverge_x*_P*_t*";

quantity = "Ux";
minValidVelocity = 0.5;

finalFraction = 0.50;          % final 50% used as reference mean
thresholds = [2.0 1.0 0.5];    % percent convergence criteria

exportFigures = false;
EXPORT_DIR = fullfile(DATA_DIR,"MeanConvergence_Output");

if exportFigures && ~exist(EXPORT_DIR,"dir")
    mkdir(EXPORT_DIR);
end

%% ================= LOAD FILES =================

files = dir(fullfile(DATA_DIR,FILE_PATTERN));

if isempty(files)
    error("No files found matching pattern: %s",FILE_PATTERN);
end

runs = struct([]);

for k = 1:numel(files)

    filepath = fullfile(files(k).folder,files(k).name);

    runs(k).file = files(k).name;
    runs(k).x = parseX(files(k).name);
    runs(k).power = parsePower(files(k).name);
    runs(k).test = parseTest(files(k).name);
    runs(k).data = readProbeFile(filepath);
    runs(k).fs = estimateSamplingFrequency(runs(k).data.t);

end

[~,order] = sortrows([[runs.x]' [runs.power]' [runs.test]']);
runs = runs(order);

fprintf("Loaded %d files.\n",numel(runs));

%% ================= ANALYZE FILES =================

results = [];
summary = table();

for k = 1:numel(runs)

    r = analyzeRun(runs(k),quantity,minValidVelocity,finalFraction,thresholds);

    if k == 1
        results = r;
    else
        results(k) = r;
    end

    newRow = table( ...
        string(r.file), ...
        r.x, ...
        r.power, ...
        r.test, ...
        r.fs, ...
        r.referenceMean, ...
        r.referenceStd, ...
        r.fluctuationPercent, ...
        r.conv2, ...
        r.conv1, ...
        r.conv05, ...
        "VariableNames", ...
        ["File", ...
         "x_mm", ...
         "Power_percent", ...
         "Test", ...
         "SamplingFrequency_Hz", ...
         "ReferenceMean_mps", ...
         "ReferenceStd_mps", ...
         "Fluctuation_percent", ...
         "Convergence_2percent_s", ...
         "Convergence_1percent_s", ...
         "Convergence_0p5percent_s"]);

    summary = [summary; newRow];

end

disp(summary)

if exportFigures
    writetable(summary,fullfile(EXPORT_DIR,"MeanConvergence_Summary.csv"));
end

%% ================= ALL RAW Ux TESTS ON ONE GRAPH =================

figure("Name","All Raw Ux Tests", ...
    "Color","w", ...
    "Position",[100 100 1200 650]);

hold on

for k = 1:numel(results)

    r = results(k);

    plot(r.t,r.y, ...
        "LineWidth",1.2, ...
        "DisplayName",sprintf("P=%g%% | Test %g",r.power,r.test));

end

hold off

grid on
box on

xlabel("Time (s)")
ylabel("U_x (m/s)")
title("Raw U_x Velocity for All Tests")

legend("Location","bestoutside")

%% ================= INDIVIDUAL TESTS VIEWER =================

labels = strings(numel(results),1);

for k = 1:numel(results)
    labels(k) = sprintf("P=%g%% | Test %g | %s", ...
        results(k).power,results(k).test,results(k).file);
end

fig1 = uifigure( ...
    "Name","Individual Mean Convergence Viewer", ...
    "Position",[100 100 1200 800]);

uilabel(fig1, ...
    "Text","Select test:", ...
    "Position",[40 755 150 25], ...
    "FontWeight","bold");

dropdown1 = uidropdown(fig1, ...
    "Items",cellstr(labels), ...
    "Value",char(labels(1)), ...
    "Position",[40 720 600 30]);

ax1 = uiaxes(fig1,"Position",[80 500 1050 180]);
ax2 = uiaxes(fig1,"Position",[80 280 1050 180]);
ax3 = uiaxes(fig1,"Position",[80 60 1050 180]);

dropdown1.ValueChangedFcn = @(src,event)plotIndividualDropdown(src,results,labels,ax1,ax2,ax3);

plotIndividualDropdown(dropdown1,results,labels,ax1,ax2,ax3);

%% ================= GROUPED BY POWER VIEWER =================

powers = unique([results.power]);

powerLabels = strings(numel(powers),1);

for k = 1:numel(powers)
    powerLabels(k) = sprintf("P=%g%%",powers(k));
end

fig2 = uifigure( ...
    "Name","Grouped Power Mean Convergence Viewer", ...
    "Position",[150 150 1200 750]);

uilabel(fig2, ...
    "Text","Select power:", ...
    "Position",[40 705 150 25], ...
    "FontWeight","bold");

dropdown2 = uidropdown(fig2, ...
    "Items",cellstr(powerLabels), ...
    "Value",char(powerLabels(1)), ...
    "Position",[40 670 200 30]);

axG1 = uiaxes(fig2,"Position",[80 390 1050 230]);
axG2 = uiaxes(fig2,"Position",[80 80 1050 230]);

dropdown2.ValueChangedFcn = @(src,event)plotGroupedDropdown(src,results,powers,powerLabels,axG1,axG2);

plotGroupedDropdown(dropdown2,results,powers,powerLabels,axG1,axG2);

%% ================= DROPDOWN VIEWER: SAME POWER COMPARISON =================

fig3 = uifigure( ...
    "Name","Same Power Test Comparison", ...
    "Position",[200 200 1200 750]);

uilabel(fig3, ...
    "Text","Select power:", ...
    "Position",[40 705 150 25], ...
    "FontWeight","bold");

dropdown3 = uidropdown(fig3, ...
    "Items",cellstr(powerLabels), ...
    "Value",char(powerLabels(1)), ...
    "Position",[40 670 200 30]);

axC1 = uiaxes(fig3,"Position",[80 390 1050 230]);
axC2 = uiaxes(fig3,"Position",[80 80 1050 230]);

dropdown3.ValueChangedFcn = @(src,event)plotSamePowerComparison(src,results,powers,powerLabels,axC1,axC2);

plotSamePowerComparison(dropdown3,results,powers,powerLabels,axC1,axC2);

%% ================================================================
%% LOCAL FUNCTIONS
%% ================================================================

function r = analyzeRun(run,quantity,minValidVelocity,finalFraction,thresholds)

T = run.data;

t = T.t;
y = T.(quantity);

good = ~isnan(t) & ~isnan(y);

if ismember(string(quantity),["Vmag","Ux","Uy","Uz"])
    good = good & y > minValidVelocity;
end

t = t(good);
y = y(good);

if numel(t) < 20
    error("Not enough valid samples in %s",run.file);
end

t = t - t(1);

N = numel(y);

referenceStartIdx = max(1,round((1-finalFraction)*N));
yReference = y(referenceStartIdx:end);

referenceMean = mean(yReference,"omitnan");
referenceStd = std(yReference,"omitnan");
fluctuationPercent = 100 * referenceStd / abs(referenceMean);

runningMean = cumsum(y) ./ (1:N)';
errorPercent = 100 * abs(runningMean - referenceMean) / abs(referenceMean);

convTimes = nan(size(thresholds));

for j = 1:numel(thresholds)

    idx = findConvergenceIndex(errorPercent,thresholds(j));

    if ~isempty(idx)
        convTimes(j) = t(idx);
    end

end

r = struct();

r.file = run.file;
r.x = run.x;
r.power = run.power;
r.test = run.test;
r.fs = run.fs;

r.t = t;
r.y = y;

r.referenceStartIdx = referenceStartIdx;
r.referenceMean = referenceMean;
r.referenceStd = referenceStd;
r.fluctuationPercent = fluctuationPercent;

r.runningMean = runningMean;
r.errorPercent = errorPercent;

r.conv2 = convTimes(thresholds == 2.0);
r.conv1 = convTimes(thresholds == 1.0);
r.conv05 = convTimes(thresholds == 0.5);

r.band2Upper = referenceMean * 1.02;
r.band2Lower = referenceMean * 0.98;

r.band1Upper = referenceMean * 1.01;
r.band1Lower = referenceMean * 0.99;

r.band05Upper = referenceMean * 1.005;
r.band05Lower = referenceMean * 0.995;

end

function idxConv = findConvergenceIndex(errorSignal,threshold)

idxConv = [];

for k = 1:length(errorSignal)

    if all(errorSignal(k:end) <= threshold)
        idxConv = k;
        return
    end

end

end

function T = readProbeFile(filepath)

lines = readlines(filepath);

rawNames = split(lines(5),sprintf("\t"));
names = cleanNames(rawNames);

opts = delimitedTextImportOptions( ...
    "Delimiter","\t", ...
    "NumVariables",numel(names), ...
    "DataLines",[6 Inf], ...
    "VariableNames",names, ...
    "VariableTypes",repmat("string",1,numel(names)));

raw = readtable(filepath,opts);

T = table();

timeStrings = raw.(names(1));

timeVals = datetime(timeStrings, ...
    "InputFormat","dd/MM/uuuu  HH:mm:ss,SSSSSS");

T.t = seconds(timeVals - timeVals(1));

for c = 2:width(raw)

    s = raw.(names(c));
    s = replace(s,",",".");

    T.(names(c)) = str2double(s);

end

end

function names = cleanNames(rawNames)

names = strings(size(rawNames));

for i = 1:numel(rawNames)

    s = string(rawNames(i));

    s = regexprep(s,"\s*\(.*?\)","");
    s = regexprep(s,"[^A-Za-z0-9]","");

    switch lower(s)

        case "time"
            s = "timeRaw";

        case "p15trel"
            s = "P1";

        case "p25trel"
            s = "P2";

        case "p35trel"
            s = "P3";

        case "p45trel"
            s = "P4";

        case "p55trel"
            s = "P5";

    end

    names(i) = matlab.lang.makeValidName(s);

end

names = matlab.lang.makeUniqueStrings(names);

end

function fs = estimateSamplingFrequency(t)

dt = median(diff(t),"omitnan");

if isempty(dt) || isnan(dt) || dt <= 0
    fs = NaN;
else
    fs = 1/dt;
end

end

function x = parseX(filename)

token = regexp(filename,"x(\d+)","tokens","once");

if isempty(token)
    x = NaN;
else
    x = str2double(token{1});
end

end

function p = parsePower(filename)

token = regexp(filename,"P(\d+)","tokens","once");

if isempty(token)
    token = regexp(filename,"p(\d+)","tokens","once");
end

if isempty(token)
    p = NaN;
else
    p = str2double(token{1});
end

end

function testNum = parseTest(filename)

token = regexp(filename,"t(\d+)","tokens","once");

if isempty(token)
    testNum = NaN;
else
    testNum = str2double(token{1});
end

end

function plotIndividualDropdown(src,results,labels,ax1,ax2,ax3)

idx = find(labels == string(src.Value),1);
r = results(idx);

cla(ax1)
plot(ax1,r.t,r.y,"LineWidth",1.0)
hold(ax1,"on")
yline(ax1,r.referenceMean,"k--","Reference mean","LineWidth",1.2)
yline(ax1,r.band1Upper,"k:","+1%","LineWidth",1.0)
yline(ax1,r.band1Lower,"k:","-1%","LineWidth",1.0)
hold(ax1,"off")

grid(ax1,"on")
box(ax1,"on")
xlabel(ax1,"Time (s)")
ylabel(ax1,"U_x (m/s)")
title(ax1,sprintf("Raw Velocity | P=%g%% | Test %g",r.power,r.test))

cla(ax2)
plot(ax2,r.t,r.runningMean,"LineWidth",1.5)
hold(ax2,"on")
yline(ax2,r.referenceMean,"k--","Reference mean","LineWidth",1.2)

yline(ax2,r.band2Upper,":","+2%","LineWidth",1.0)
yline(ax2,r.band2Lower,":","-2%","LineWidth",1.0)
yline(ax2,r.band1Upper,":","+1%","LineWidth",1.0)
yline(ax2,r.band1Lower,":","-1%","LineWidth",1.0)
yline(ax2,r.band05Upper,":","+0.5%","LineWidth",1.0)
yline(ax2,r.band05Lower,":","-0.5%","LineWidth",1.0)

if ~isnan(r.conv2)
    xline(ax2,r.conv2,"--","2%","LineWidth",1.2)
end

if ~isnan(r.conv1)
    xline(ax2,r.conv1,"--","1%","LineWidth",1.2)
end

if ~isnan(r.conv05)
    xline(ax2,r.conv05,"--","0.5%","LineWidth",1.2)
end

hold(ax2,"off")

grid(ax2,"on")
box(ax2,"on")
xlabel(ax2,"Time (s)")
ylabel(ax2,"Running mean U_x (m/s)")
title(ax2,sprintf("Mean Convergence | 2%% = %.2f s | 1%% = %.2f s | 0.5%% = %.2f s", ...
    r.conv2,r.conv1,r.conv05))

cla(ax3)
plot(ax3,r.t,r.errorPercent,"LineWidth",1.5)
hold(ax3,"on")
yline(ax3,2,"--","2%","LineWidth",1.0)
yline(ax3,1,"--","1%","LineWidth",1.0)
yline(ax3,0.5,"--","0.5%","LineWidth",1.0)

if ~isnan(r.conv2)
    xline(ax3,r.conv2,"--","2%","LineWidth",1.2)
end

if ~isnan(r.conv1)
    xline(ax3,r.conv1,"--","1%","LineWidth",1.2)
end

if ~isnan(r.conv05)
    xline(ax3,r.conv05,"--","0.5%","LineWidth",1.2)
end

hold(ax3,"off")

grid(ax3,"on")
box(ax3,"on")
xlabel(ax3,"Time (s)")
ylabel(ax3,"Running mean error (%)")
title(ax3,"Running Mean Error Relative to Reference Mean")

end

function plotGroupedDropdown(src,results,powers,powerLabels,ax1,ax2)

idxPower = find(powerLabels == string(src.Value),1);
p = powers(idxPower);

idx = find([results.power] == p);
rset = results(idx);

minEndTime = min(arrayfun(@(r) r.t(end),rset));
tCommon = linspace(0,minEndTime,1000)';

errorMatrix = nan(numel(tCommon),numel(rset));

for j = 1:numel(rset)
    errorMatrix(:,j) = interp1(rset(j).t,rset(j).errorPercent,tCommon,"linear");
end

avgError = mean(errorMatrix,2,"omitnan");

conv2_all = [rset.conv2];
conv1_all = [rset.conv1];
conv05_all = [rset.conv05];

valid2 = ~isnan(conv2_all);
valid1 = ~isnan(conv1_all);
valid05 = ~isnan(conv05_all);

avgConv2 = mean(conv2_all(valid2));
avgConv1 = mean(conv1_all(valid1));
avgConv05 = mean(conv05_all(valid05));

n2 = sum(valid2);
n1 = sum(valid1);
n05 = sum(valid05);

nTotal = numel(rset);

cla(ax1)

plot(ax1,tCommon,errorMatrix,"LineWidth",1.0)
hold(ax1,"on")
plot(ax1,tCommon,avgError,"k","LineWidth",2.5)

yline(ax1,2,"--","2%","LineWidth",1.0)
yline(ax1,1,"--","1%","LineWidth",1.0)
yline(ax1,0.5,"--","0.5%","LineWidth",1.0)

if ~isnan(avgConv2)
    xline(ax1,avgConv2,"--","Mean 2%","LineWidth",1.2)
end

if ~isnan(avgConv1)
    xline(ax1,avgConv1,"--","Mean 1%","LineWidth",1.2)
end

if ~isnan(avgConv05)
    xline(ax1,avgConv05,"--","Mean 0.5%","LineWidth",1.2)
end

hold(ax1,"off")

grid(ax1,"on")
box(ax1,"on")
xlabel(ax1,"Time (s)")
ylabel(ax1,"Running mean error (%)")
title(ax1,sprintf("Grouped Mean Convergence | P=%g%% | Mean: 2%% = %.2f s | 1%% = %.2f s | 0.5%% = %.2f s", ...
    p,avgConv2,avgConv1,avgConv05))

legend(ax1,["Test 1","Test 2","Test 3","Average"],"Location","northeast")

cla(ax2)

axis(ax2,[0 1 0 1])
axis(ax2,'off')

text(ax2,0.05,0.80,...
    sprintf('Power Level = %.0f%%',p),...
    'FontSize',16,...
    'FontWeight','bold')

text(ax2,0.05,0.60,...
    sprintf('Average 2%% Convergence Time: %.2f s (%d/%d reached)',avgConv2,n2,nTotal),...
    'FontSize',14)

text(ax2,0.05,0.42,...
    sprintf('Average 1%% Convergence Time: %.2f s (%d/%d reached)',avgConv1,n1,nTotal),...
    'FontSize',14)

text(ax2,0.05,0.24,...
    sprintf('Average 0.5%% Convergence Time: %.2f s (%d/%d reached)',avgConv05,n05,nTotal),...
    'FontSize',14)

text(ax2,0.05,0.08,...
    sprintf('Based on %d repeated tests',numel(rset)),...
    'FontSize',12)

end

function plotSamePowerComparison(src,results,powers,powerLabels,ax1,ax2)

idxPower = find(powerLabels == string(src.Value),1);
p = powers(idxPower);

idx = find([results.power] == p);
rset = results(idx);

cla(ax1)
hold(ax1,"on")

for j = 1:numel(rset)
    plot(ax1,rset(j).t,rset(j).runningMean, ...
        "LineWidth",1.5, ...
        "DisplayName",sprintf("Test %g",rset(j).test));
end

avgRef = mean([rset.referenceMean],"omitnan");

yline(ax1,avgRef,"k--","Average reference mean","LineWidth",1.3)

hold(ax1,"off")

grid(ax1,"on")
box(ax1,"on")
xlabel(ax1,"Time (s)")
ylabel(ax1,"Running mean U (m/s)")
title(ax1,sprintf("Running Mean Comparison at Same Power | P=%g%%",p))
legend(ax1,"Location","best")

% Force same x scale across all same-power tests
maxTime = max(arrayfun(@(r) r.t(end),rset));
xlim(ax1,[0 maxTime])

cla(ax2)
hold(ax2,"on")

for j = 1:numel(rset)
    plot(ax2,rset(j).t,rset(j).errorPercent, ...
        "LineWidth",1.5, ...
        "DisplayName",sprintf("Test %g",rset(j).test));
end

yline(ax2,2,"--","2%","LineWidth",1.0)
yline(ax2,1,"--","1%","LineWidth",1.0)
yline(ax2,0.5,"--","0.5%","LineWidth",1.0)

hold(ax2,"off")

grid(ax2,"on")
box(ax2,"on")
xlabel(ax2,"Time (s)")
ylabel(ax2,"Running mean error (%)")
title(ax2,sprintf("Running Mean Error Comparison at Same Power | P=%g%%",p))
legend(ax2,"Location","best")

xlim(ax2,[0 maxTime])

end