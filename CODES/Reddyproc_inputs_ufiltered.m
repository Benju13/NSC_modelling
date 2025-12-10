%%%% Last modified: 2024-01-23 %%%%% 
%%% Purpose: Process AmeriFlux data with ustar filtering and prepare for ReddyProc
% No gapfilling applied - uses real data only
% This includes ustar filtering 
% If you remove neg respiration then NEE becomes 150 and then almost equal
% Input: AmeriFlux CSV files
% Output: ReddyProc-formatted CSV files

%%
clear;
clc;

%% Extract annaul data
global yyear
sitename = 'CRK';
yyear = 2024;

if strcmp(sitename, "CST")
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\CST\AMF_US-Cst_BASE-BADM_1-5\US-Cst_HR_201301010000_201712312330_v4.csv");
elseif strcmp(sitename, "GA")
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\GA\AMF_US-LL1_BASE_HH_2-5.csv");
elseif strcmp(sitename, "CRK") & yyear == 2022
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\CRK\US-CRK_HH_202201010000_202301010000.csv");
    elseif strcmp(sitename, "CRK") & yyear == 2023
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\CRK\US-CRK_HH_202301010000_202401010000.csv");
    elseif strcmp(sitename, "CRK") & yyear == 2024
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\CRK\US-CRK_HH_202401010000_202501010000.csv");
    elseif strcmp(sitename, "CRK") & yyear == 2025
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\CRK\US-CRK_HH_202501010000_202507010000.csv");

elseif strcmp(sitename, "NC2")
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\NC2\AMF_US-NC2_BASE_HH_11-5.csv");
end
extract_year(data,yyear,sitename)
close all; 

%% Inputs

fclose('all');
sitename = 'CRK';
if strcmp(sitename, 'CST')
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\CST\AMF_US-Cst_BASE-BADM_1-5\US-Cst_HR_201301010000_201712312330_v4.csv");
elseif strcmp(sitename, 'CRK') && yyear == 2023
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\CRK\US-CRK_HH_202301010000_202401010000.csv");
elseif strcmp(sitename, 'CRK') && yyear == 2022
   % data = readtable("C:\Benju\CRK_processing\Output_data\Final data product for Ameriflux submission\US-CRK_HH_202201010000_202301010000.csv");
data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\CRK\US-CRK_HH_202201010000_202301010000.csv");
elseif strcmp(sitename, 'CRK') && yyear == 2024
data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\CRK\US-CRK_HH_202401010000_202501010000.csv");
elseif strcmp(sitename, 'GA')
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\GA\AMF_US-LL1_BASE_HH_2-5.csv");
    elseif strcmp(sitename, 'NC2')
    data = readtable("C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ameriflux\NC2\AMF_US-NC2_BASE_HH_11-5.csv");
end

%%
%%% Datetime change %%%
x=num2str(data.TIMESTAMP_END);
data.dt= datetime(x,'InputFormat','yyyyMMddHHmm','Format','yyyy/MM/dd HH:mm');

ts=(datenum(num2str(data.TIMESTAMP_START),'yyyymmddHHMM'));
doy = day(datetime(ts,'ConvertFrom','datenum'),'dayofyear');
Year = year(datetime(ts, 'ConvertFrom','datenum'),'gregorian');
hours=hour(datetime(ts,'ConvertFrom','datenum'))+minute(datetime(ts,'ConvertFrom','datenum'))/60;
data.doy = doy; 
data.hours = hours;
data.Year = Year;
data.Site = repmat(sitename,size(data,1),1);

%Filter data by year 
index = (data.Year == yyear);
data = data(index,:);

% Add month snd season : 
months = month(data.dt);
spring_ind = (months == 3) | (months == 4) | (months == 5);
summer_ind = (months == 6) | (months == 7) | (months == 8);
fall_ind = (months == 9) | (months == 10) | (months == 11);
winter_ind = (months == 12) | (months == 1) | (months == 2);

data.season = cell(size(data,1),1);
data.season(spring_ind) = {'spring'};
data.season(summer_ind) = {'summer'};
data.season(fall_ind) = {'fall'};
data.season(winter_ind) = {'winter'};
data.season = char(data.season);

%% Ustar data
directory ='C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\Ustar\';
filename = strcat('Ustar_', sitename, '.xlsx');
fullFilePath = fullfile(directory, filename);
ustar_data = readtable(fullFilePath);
% Read ustar data 
ustar_data.Site = categorical(ustar_data.Site);  % Convert to categorical array
ustar_data.Site = char(ustar_data.Site);  % Convert to character array
ustar_data.season = categorical(ustar_data.season);  % Convert to categorical array
ustar_data.season = char(ustar_data.season);  % Convert to character array

data = outerjoin(data, ustar_data, 'Keys', {'Site', 'Year', 'season'}, 'Type', 'left', 'MergeKeys', true);

% Sort the merged data based on the timestamp column (assuming it's named 'timestamp')
data = sortrows(data, 'TIMESTAMP_START');
data = unique(data);

%Change to Nan for plotting 
varNames = data.Properties.VariableNames;
for i = 1:numel(varNames)
    variable = data.(varNames{i});
    if isnumeric(variable)
        variable(variable == -9999) = NaN;
    end
    data.(varNames{i}) = variable;
end

% Add storage fluxes in CST data 
if strcmp(sitename, 'CST') 
    for i = 1:size(data.FC)
    if ~isnan(data.FC(i))
    data.NEE(i) = data.FC(i) + data.SC(i);
    else data.NEE(i) = NaN;
end
    end 
end

%% NC2 data 
if strcmp(sitename, 'NC2')  
if yyear == 2011 || yyear ==2012 || yyear == 2013 || yyear == 2014 || yyear == 2015 || yyear == 2016 || yyear == 2017 || yyear == 2018|| yyear == 2019
    data.NEE_1_1_1 = data.NEE_PI_1_1_1;
elseif yyear == 2008 || yyear == 2009 || yyear == 2010
    data.NEE_1_1_1 = data.NEE_PI_F_1_1_1;
end
end
% Ustar filtering 
if strcmp(sitename,'CRK')
    filter_condition= data.USTAR_1_1_1 < data.ustar_thres ; 
data(filter_condition, 'NEE_1_1_1') = {NaN};
elseif strcmp(sitename,'CST')
filter_condition= data.USTAR < data.ustar_thres ; 
data(filter_condition, 'NEE') = {NaN};
elseif strcmp(sitename,'GA')
filter_condition= data.USTAR < data.ustar_thres ; 
data(filter_condition, 'NEE_PI_F') = {NaN};
elseif strcmp(sitename,'NC2')
filter_condition= data.USTAR_1_1_1 < data.ustar_thres ; 
data(filter_condition, 'NEE_1_1_1') = {NaN};
end

%% GA
%Apply quality control for SWIN Look at April 1 to April 13
if strcmp(sitename, 'GA') && yyear == 2019 
start_date = datetime(2019, 4, 2);
end_date = datetime(2019, 4, 18);
rows_to_delete = data.dt >= start_date & data.dt <= end_date;
data.SW_IN_PI_F (rows_to_delete, :) = NaN;
end

%% CST
if strcmp(sitename, 'CST') 
    data.Reco_org(data.Reco_org > 20 | data.Reco_org < -0.1) = NaN;
    data.NEE(data.time == 0) = data.Reco_org(data.time == 0);
end

%% Averaginf soil data 
if strcmp(sitename, 'CST') 
data.Tsoil_avg = (data.TS_1_1_1);
elseif strcmp(sitename, 'CRK') 
data.Tsoil_avg = (data.TS_1_1_1+data.TS_1_2_1)/2;
elseif strcmp(sitename, 'GA')
data.Tsoil_avg = (data.TS_1_1_1+data.TS_1_2_1)/2;
elseif strcmp(sitename, 'NC2')
data.Tsoil_avg = (data.TS_1_1_1+data.TS_1_2_1)/2;
end

%TESTSSS
for i = 1:size(data.PPFD_IN_1_1_1)
if data.PPFD_IN_1_1_1 (i) < 5
    data.SW_IN_1_1_1(i) = 0;
end
end
%% FInal data 
if strcmp(sitename, 'GA')
    Year = data.Year;
DoY = data.doy;
Hour = data.hours;
NEE = data.NEE_PI_F;
H = data.H_PI_F;
LE = data.LE_PI_F; 
Rg = data.SW_IN_PI_F;
Tair = data.TA;
Tsoil=data.Tsoil_avg;
rH = data.RH;
VPD = data.VPD_PI;
Ustar =data.USTAR;
PPFD = data.PPFD_IN_PI_F;
rpdata = table(Year,DoY,Hour,NEE,LE,H,Rg,Tair,Tsoil,rH,VPD,Ustar,PPFD);

elseif strcmp(sitename, 'CST')
Year = data.Year;
DoY = data.doy;
Hour = data.hours;
NEE = data.NEE;
H = data.H;
LE = data.LE; 
Rg = data.SW_IN;
Tair = data.Ta;
Tsoil=data.Tsoil_avg;
rH = data.RH;
VPD = data.VPD;
Ustar =data.USTAR;
PPFD = data.PPFD_IN;
rpdata = table(Year,DoY,Hour,NEE,LE,H,Rg,Tair,Tsoil,rH,VPD,Ustar,PPFD);

elseif strcmp(sitename, 'CRK')
Year = data.Year;
DoY = data.doy;
Hour = data.hours;
NEE = data.NEE_1_1_1;
H = data.H_1_1_1;
LE = data.LE_1_1_1; 
Rg = data.SW_IN_1_1_1;
Tair = data.TA_1_1_1;
Tsoil=data.Tsoil_avg;
rH = data.RH_1_1_1;
VPD = data.VPD_1_1_1;
Ustar =data.USTAR_1_1_1;
PPFD = data.PPFD_IN_1_1_1;
%Added by Benju on 2024-07-15
SWC_5cm = data.SWC_1_1_1;
SWC_20cm = data.SWC_1_2_1;
QC_flag = data.FC_SSITC_TEST_1_1_1;

%% Removal of extreme negative nightime values based on precipitation and storm event in CRK
if strcmp(sitename, 'CRK') && yyear == 2022
% Define specific problematic days
%Doy 54 to doy 57 correlate with rain event 
%34 with 
%67 also lot of rain 
    problematic_days = [34, 54, 55, 56, 57, 58,67];
    is_problematic_day = ismember(DoY, problematic_days);
    NEE(is_problematic_day) = NaN;
    data.NEE_1_1_1(is_problematic_day) = NaN;
    data.FC_1_1_1(is_problematic_day) = NaN;
    num_filtered = sum(is_problematic_day);
    disp(['Set ', num2str(num_filtered), ' NEE values to NaN for DoY ', ...
          num2str(problematic_days), ' at CRK 2022.']);
end
rpdata = table(Year,DoY,Hour,NEE,LE,H,Rg,Tair,Tsoil,rH,VPD,Ustar,PPFD,SWC_5cm,SWC_20cm,QC_flag);

elseif strcmp(sitename, 'NC2')
Year = data.Year;
DoY = data.doy;
Hour = data.hours;
NEE = data.NEE_1_1_1;
H = data.H_1_1_1;
LE = data.LE_1_1_1; 
Rg = data.SW_IN_1_1_1;
Tair = data.TA_1_1_1;
Tsoil=data.Tsoil_avg;
rH = data.RH_1_1_1;
VPD = data.VPD_PI_1_1_1;
Ustar =data.USTAR_1_1_1;
PPFD = data.PPFD_IN_1_1_1;
rpdata = table(Year,DoY,Hour,NEE,LE,H,Rg,Tair,Tsoil,rH,VPD,Ustar,PPFD);
end

directory = ['C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\Input Data\ReddyProc\Seasonal\',sitename];
filename = strcat('rp',sitename, '_', num2str(yyear), '.csv');
fullFilePath = fullfile(directory, filename);
writetable(rpdata, fullFilePath);
disp(['Code executed for site: ', sitename, ', year: ', num2str(yyear)]);


%% Plotting nightime NEE
is_night = (PPFD < 5) | (Rg <= 20);
nighttime_nee = NEE(is_night);
nighttime_doy = DoY(is_night);

% Plot
figure;
scatter(nighttime_doy, nighttime_nee, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Day of Year');
ylabel('Nighttime NEE (μmol CO₂ m⁻² s⁻¹)');
title(['Nighttime NEE - ', sitename, ' ', num2str(yyear)]);
grid on;
xlim([1 366]);
yline(0, 'k--', 'LineWidth', 1);

%% More figure 
%Diagnostic plot: Nighttime NEE vs u* by temperature bins
is_night = (PPFD < 5) | (Rg <= 20);
nee_night   = NEE(is_night);
ustar_night = Ustar(is_night);
tair_night  = Tair(is_night);

% Remove NaNs
valid = ~isnan(nee_night) & ~isnan(ustar_night) & ~isnan(tair_night);
nee_night   = nee_night(valid);
ustar_night = ustar_night(valid);
tair_night  = tair_night(valid);

% Define temperature bins (quartiles here; change nbins if needed)
nbins = 4;
edges = quantile(tair_night, linspace(0,1,nbins+1));
[~,~,binID] = histcounts(tair_night, edges);

% Colors for bins
cmap = lines(nbins);

figure; hold on;
for b = 1:nbins
    idx = binID == b;
    if sum(idx) > 10
        scatter(ustar_night(idx), nee_night(idx), 12, cmap(b,:), 'filled', 'MarkerFaceAlpha',0.3);
        % Bin-averaged trend
        [ustar_sorted, sortIdx] = sort(ustar_night(idx));
        nee_sorted = nee_night(idx);
        % smooth moving average
        win = max(5, round(numel(sortIdx)/20));
        nee_smooth = movmean(nee_sorted(sortIdx), win);
        plot(ustar_sorted, nee_smooth, 'Color', cmap(b,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('T %.1f–%.1f °C', edges(b), edges(b+1)));
    end
end
xlabel('u* (m s^{-1})');
ylabel('Nighttime NEE (\mumol CO_{2} m^{-2} s^{-1})');
title(sprintf('%s %d Nighttime NEE vs u* by Tair bins', sitename, yyear));
yline(0,'k--'); grid on;
legend('show','Location','best');

