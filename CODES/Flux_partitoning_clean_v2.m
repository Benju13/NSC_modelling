%═══════════════════════════════════════════════════════════════════════════
% REDDYPROC NIGHTTIME PARTITIONING CORRECTION - BATCH PROCESSOR
%═══════════════════════════════════════════════════════════════════════════
%
% PURPOSE:
%   This script processes gap-filled eddy covariance data from REddyProc to
%   apply nighttime partitioning corrections for GPP and Reco calculations.
%   It handles multiple sites and years automatically.
%
% WHAT IT DOES:
%   1. Reads REddyProc gap-filled output files (NEE, Reco, GPP)
%   2. Merges sunrise/sunset data to identify day/night periods
%   3. Applies nighttime corrections:
%      - Sets nighttime GPP = 0 (no photosynthesis in darkness)
%      - Uses measured NEE as Reco during nighttime (when available)
%      - Uses modeled Reco during daytime
%   4. Recalculates corrected GPP from mass balance: GPP = Reco - NEE
%   5. Outputs corrected files with "_corrected.txt" suffix
%
% INPUTS REQUIRED:
%   - REddyProc output files: [SITE]_[YEAR]_rpresult.txt
%   - Sunrise/sunset files: [SITE]_navy_[YEAR].xlsx
%
% OUTPUTS:
%   - Corrected files: [SITE]_[YEAR]_rpresult_corrected.txt
%   - Summary statistics printed to console
%
% USER CONFIGURATION:
%   Edit the following sections below:
%   - ustar: Choose "seasonal" or "reddyproc" for u* filtering method
%   - sites: List of site codes to process (e.g., {'GA', 'NC2', 'CRK'})
%   - years: Range of years to process (e.g., 2008:2024)
%   - base_dir: Main directory path for your data
%
% AUTHOR: Benju Baniya
% LAST UPDATED: November 14, 2025
%═══════════════════════════════════════════════════════════════════════════

%% ═══════════════════════════════════════════════════════════════════════
%  USER CONFIGURATION - EDIT THIS SECTION
%  ═══
clear; clc;

% Configuration
ustar = "seasonal";  % or "reddyproc"or "seasonal"

% Define sites and years to process
sites = {'GA', 'NC2', 'CRK', 'CST'};
years = 2022:2024;  % Adjust as needed

base_dir = 'C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\';
input_dir = fullfile(base_dir, 'Input Data');
output_dir = fullfile(base_dir, 'Output Data');

%% Loop through all sites and years
for s = 1:length(sites)
    sitename = sites{s};
    
    for y = 1:length(years)
        year = years(y);
        
        fprintf('\n========================================\n');
        fprintf('Processing: %s - %d\n', sitename, year);
        fprintf('========================================\n');
        
        try
            % Set paths based on ustar type
            if strcmp(ustar, "seasonal")
                reddyproc_path = fullfile(output_dir, 'Reddyproc', 'Seasonal', sitename);
            else
                reddyproc_path = fullfile(output_dir, 'Reddyproc', 'Reddyproc_ustar', sitename);
            end
            
            %% 1) Read sunrise/sunset
            sunrise_path = fullfile(input_dir, 'sunrise_sunset', sitename);
            sun_file = fullfile(sunrise_path, sprintf('%s_navy_%d.xlsx', sitename, year));
            
            if ~isfile(sun_file)
                warning('Skipping %s %d: Sunrise file not found', sitename, year);
                continue;
            end
            
            sunrise_sunset = readtable(sun_file);
            
            % Compute DoY if needed
            if ~ismember('DoY', sunrise_sunset.Properties.VariableNames) && ...
               all(ismember({'Month','Day'}, sunrise_sunset.Properties.VariableNames))
                sunrise_sunset.DoY = day(datetime(year, sunrise_sunset.Month, sunrise_sunset.Day), 'dayofyear');
            end
            
            % Convert Month to numeric if needed
            if ischar(sunrise_sunset.Month) || isstring(sunrise_sunset.Month)
                sunrise_sunset.Month = month(datetime(sunrise_sunset.Month,'InputFormat','MMM'));
            end
            
            % Parse sunrise/sunset times
            sunrise_sunset.Sunrise_new = datenum(num2str(sunrise_sunset.Sunrise), 'HHMM');
            sunrise_sunset.Sunset_new  = datenum(num2str(sunrise_sunset.Sunset), 'HHMM');
            sunrise_sunset.Year = repmat(year, height(sunrise_sunset), 1);
            
            %% 2) Read gap-filled ReddyProc output
            rp_file = fullfile(reddyproc_path, sprintf('%s_%d_rpresult.txt', sitename, year));
            
            if ~isfile(rp_file)
                warning('Skipping %s %d: ReddyProc file not found', sitename, year);
                continue;
            end
            
            data = readtable(rp_file);
            
            % Special case for GA 2018
            if strcmp(sitename,'GA') && year==2018
                data.Reco(13894:end) = -9999;
            end
            
            %% 3) Build datetime, Month, DoY
            data.dt = datetime(data.Year,1,0) + days(data.DoY) + hours(data.Hour);
            data.Month = month(data.dt);
            
            %% 4) Merge sunrise/sunset
            data = join(data, sunrise_sunset(:,{'Year','DoY','Sunrise_new','Sunset_new'}), 'Keys', {'Year','DoY'});
            
            %% 5) Clean up -9999 → NaN
            vars = data.Properties.VariableNames;
            for i = 1:numel(vars)
                if isnumeric(data.(vars{i}))
                    data.(vars{i})(data.(vars{i}) == -9999) = NaN;
                end
            end
            
            %% 6) Flag day vs. night
            n = height(data);
            sunriseHour = (data.Sunrise_new - floor(data.Sunrise_new)) * 24;
            sunsetHour  = (data.Sunset_new - floor(data.Sunset_new)) * 24;
            timeHour = hour(data.dt) + minute(data.dt)/60;
            isNight_sun = (timeHour < sunriseHour) | (timeHour > sunsetHour);
            isNight_ppfd = data.PPFD < 4;
            data.IsNight = isNight_sun | isNight_ppfd;
            data.IsDay = ~data.IsNight;
            
            %% 7-8) Build Reco_NT_f and GPP_NT_f
            data.Reco_NT_f = nan(n,1);
            data.GPP_NT_f = nan(n,1);
            
            if strcmp(ustar, "seasonal")
                % Seasonal ustar approach
                m1 = data.IsNight & ~isnan(data.NEE_orig);
                data.Reco_NT_f(m1) = data.NEE_orig(m1);
                
                m2 = data.IsNight & isnan(data.NEE_orig);
                data.Reco_NT_f(m2) = data.NEE_f(m2);
                
                m3 = data.IsDay;
                data.Reco_NT_f(m3) = data.Reco(m3);
                
                % GPP
                data.GPP_NT_f(data.IsNight) = 0;
                data.GPP_NT_f(data.IsDay) = data.Reco_NT_f(data.IsDay) - data.NEE_f(data.IsDay);
                
            else
                % ReddyProc ustar approach
                 % Multiple uncertainty estimates
                data.Reco_NT_f_u05 = nan(n,1);
                data.Reco_NT_f_u50 = nan(n,1);
                data.Reco_NT_f_u95 = nan(n,1);

                m1 = data.IsNight & ~isnan(data.NEE_uStar_orig);
                data.Reco_NT_f_u05(m1) = data.NEE_uStar_orig(m1);
                data.Reco_NT_f_u50(m1) = data.NEE_uStar_orig(m1);
                data.Reco_NT_f_u95(m1) = data.NEE_uStar_orig(m1);
                
                m2 = data.IsNight & isnan(data.NEE_uStar_orig);
                data.Reco_NT_f_u05(m2) = data.NEE_uStar_f(m2);
                data.Reco_NT_f_u50(m2) = data.NEE_uStar_f(m2);
                data.Reco_NT_f_u95(m2) = data.NEE_uStar_f(m2);
                
                m3 = data.IsDay;
                data.Reco_NT_f_u05(m3) = data.Reco_U05(m3);
                data.Reco_NT_f_u50(m3) = data.Reco_U50(m3);
                data.Reco_NT_f_u95(m3) = data.Reco_U95(m3);
                
                % GPP
                data.GPP_NT_f_u05 = nan(n,1);
                data.GPP_NT_f_u50 = nan(n,1);
                data.GPP_NT_f_u95 = nan(n,1);
                
                data.GPP_NT_f_u05(data.IsNight) = 0;
                data.GPP_NT_f_u50(data.IsNight) = 0;
                data.GPP_NT_f_u95(data.IsNight) = 0;
                
                data.GPP_NT_f_u05(data.IsDay) = data.Reco_NT_f_u05(data.IsDay) - data.NEE_uStar_f(data.IsDay);
                data.GPP_NT_f_u50(data.IsDay) = data.Reco_NT_f_u50(data.IsDay) - data.NEE_uStar_f(data.IsDay);
                data.GPP_NT_f_u95(data.IsDay) = data.Reco_NT_f_u95(data.IsDay) - data.NEE_uStar_f(data.IsDay);
            end
            
            %% 9) Summary statistics
            fprintf('\n=== NIGHTTIME CORRECTION SUMMARY ===\n');
            if strcmp(ustar, "seasonal")
                night_measured = sum(data.IsNight & ~isnan(data.NEE_orig));
            else
                night_measured = sum(data.IsNight & ~isnan(data.NEE_uStar_orig));
            end
            night_total = sum(data.IsNight);
            fprintf('Total nighttime periods: %d\n', night_total);
            fprintf('Nighttime with measured NEE: %d (%.1f%%)\n', night_measured, 100*night_measured/night_total);
            
            %% 10) Write output
            out_name = sprintf('%s_%d_rpresult_corrected.txt', sitename, year);
            out_path = fullfile(reddyproc_path, out_name);
            writetable(data, out_path, 'Delimiter', '\t', 'FileType', 'text', 'QuoteStrings', false);
            fprintf('✓ Wrote: %s\n', out_name);
            
        catch ME
            warning('Error processing %s %d: %s', sitename, year, ME.message);
            continue;
        end
    end
end

fprintf('\n========================================\n');
fprintf('BATCH PROCESSING COMPLETE\n');
fprintf('========================================\n');