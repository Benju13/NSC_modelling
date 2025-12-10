%═══════════════════════════════════════════════════════════════════════════
% SEASONAL REDDYPROC ANALYSIS - CREATING MONTHLY AND DAILY DATA
%═══════════════════════════════════════════════════════════════════════════
%
% PURPOSE:
%   Processes seasonal REddyProc data (day/night partitioning) for multiple 
%   sites and years. Creates diagnostic figures comparing DT and NT methods.
%
% OUTPUTS:
%   - Monthly and daily CSV files for each site-year
%   - Diagnostic PNG figures for respiration models
%   - DT vs NT comparison plots
%   - Summary statistics printed to console
%
% AUTHOR: Benju Baniya
% LAST UPDATED: November 16, 2025
%═══════════════════════════════════════════════════════════════════════════

%% ═══════════════════════════════════════════════════════════════════════
%  USER CONFIGURATION
%  ═══════════════════════════════════════════════════════════════════════

clear; clc;

% Define sites and years to process
sites = {'NC2', 'GA', 'CRK', 'CST'};
years = 2022:2024;

% Conversion multiplier for flux units
multiplier = 12*60*60*24*10^(-6);
secs_per_interval = 60 * 30; % 30-minute intervals

% Base directories
base_dir = 'C:\Benju\Matlab_data_play\SouthernPine_DataAnalysis\';
config = struct();
config.multiplier = multiplier;
config.base_dir = base_dir;
config.output_dir = fullfile(base_dir, 'Output Data');
config.input_dir = fullfile(base_dir, 'Input Data');
config.fig_dir = fullfile(config.output_dir, 'Figures', 'Model Comparision', 'Seasonal');

% Create output directory if it doesn't exist
if ~exist(config.fig_dir, 'dir')
    mkdir(config.fig_dir);
end

%% ═══════════════════════════════════════════════════════════════════════
%  BATCH PROCESSING LOOP
%  ═══════════════════════════════════════════════════════════════════════

fprintf('\n╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║  SEASONAL REDDYPROC ANALYSIS - BATCH PROCESSOR                ║\n');
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');

% Track statistics
total_processed = 0;
total_skipped = 0;
total_errors = 0;

% Process each site
for s = 1:length(sites)
    sitename = sites{s};
    
    fprintf('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
    fprintf('PROCESSING SITE: %s\n', sitename);
    fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n');
    
    for y = 1:length(years)
        yyear = years(y);
        
        fprintf('  Year %d (%d/%d)... ', yyear, y, length(years));
        
        try
            %% STEP 1: READ SUNRISE/SUNSET DATA
            directory = [config.input_dir, '\sunrise_sunset\', sitename];
            filename = sprintf('%s_navy_%d.xlsx', sitename, yyear);
            fullFilePath = fullfile(directory, filename);
            
            if ~isfile(fullFilePath)
                fprintf('⚠ SKIPPED (sunrise file not found)\n');
                total_skipped = total_skipped + 1;
                continue;
            end
            
            sunrise_sunset = readtable(fullFilePath);
            sunrise_sunset.Month = char(sunrise_sunset.Month);
            sunrise_sunset.Sunrise = num2str(sunrise_sunset.Sunrise);
            sunrise_sunset.Sunset = num2str(sunrise_sunset.Sunset);
            sunrise_sunset.Sunrise_new = datenum(sunrise_sunset.Sunrise, 'HHMM');
            sunrise_sunset.Sunset_new = datenum(sunrise_sunset.Sunset, 'HHMM');
            
            %% STEP 2: READ SEASONAL REDDYPROC DATA
            if strcmp(sitename, 'NC2')
                directory = [config.output_dir, '\Reddyproc\Reddyproc_ustar\', sitename];
            else
                directory = [config.output_dir, '\Reddyproc\Seasonal\', sitename];
            end
            
            filename = sprintf('%s_%d_rpresult_corrected.txt', sitename, yyear);
            fullFilePath = fullfile(directory, filename);
            
            if ~isfile(fullFilePath)
                fprintf('⚠ SKIPPED (REddyProc file not found)\n');
                total_skipped = total_skipped + 1;
                continue;
            end
            
            data = readtable(fullFilePath);
            
            %% STEP 3: CHECK FOR DUPLICATE TIMESTAMPS
            temp_dt = datetime(data.Year, 1, 0) + days(data.DoY) + hours(data.Hour);
            [~, unique_idx] = unique(temp_dt, 'stable');
            if length(unique_idx) < height(data)
                fprintf('⚠ Warning: %d duplicate timestamps found and removed. ', height(data) - length(unique_idx));
                data = data(unique_idx, :);
            end
            
            data.dt = datetime(data.Year, 1, 0) + days(data.DoY) + hours(data.Hour);
            
            %% STEP 4: SPECIAL CASES
            % GA 2018 fix
            if strcmp(sitename, 'GA') && yyear == 2018
                data.Reco(13894:end) = -9999;
            end
            
            % CRK 2023 fix - delete OND Reco
            if strcmp(sitename, 'CRK') && yyear == 2023
                data.Reco(12827:end) = -9999;
            end
            
            %% STEP 5: PROCESS DATETIME AND MONTHS
            months = {'Jan','Feb','March','April','May','June','July','August','September','Oct','Nov','Dec'};
            month_nums = month(data.dt);
            data.Month = transpose(months(month_nums));
            data.Month = char(data.Month);
            
            %% STEP 6: MERGE SUNRISE/SUNSET DATA
            % Convert Month to char for both tables
            sunrise_sunset.Month = char(sunrise_sunset.Month);
            data.Month = char(data.Month);
            
            % Simple outerjoin like original code
            leftKeys = {'Year', 'DoY', 'Month'};
            rightKeys = {'Year', 'DoY', 'Month'};
            data = outerjoin(data, sunrise_sunset, 'LeftKeys', leftKeys, 'RightKeys', rightKeys);
            
            % Process sunrise/sunset times
            data.Sunrise = num2str(data.Sunrise);
            data.Sunset = num2str(data.Sunset);
            data.Sunrise = datenum(data.Sunrise, 'HHMM');
            data.Sunset = datenum(data.Sunset, 'HHMM');
            
            %% STEP 7: CREATE DAY/NIGHT FLAGS
            hourss = data.dt.Hour;
            minn = data.dt.Minute;
            time_to_check = datenum(hourss+":"+minn, 'HH:MM');
            data.time = zeros(height(data), 1);
            data.time(time_to_check >= data.Sunrise & time_to_check <= data.Sunset) = 1;
            
            % Create Reco_org (measured nighttime NEE)
            data.Reco_org = NaN(height(data), 1);
            
            night_idx = data.time == 0;
            
                height(data), length(night_idx), sum(night_idx);
            
            if ~ismember('NEE', data.Properties.VariableNames)
                error('NEE column does not exist in data table');
            end
            
            nighttime_nee = data.NEE(night_idx);
            
            data.Reco_org(night_idx) = nighttime_nee;
            
            %% STEP 8: REPLACE -9999 WITH NaN
            varNames = data.Properties.VariableNames;
            for i = 1:numel(varNames)
                variable = data.(varNames{i});
                if isnumeric(variable)
                    variable(variable == -9999) = NaN;
                end
                data.(varNames{i}) = variable;
            end
            
            %% STEP 9: CREATE NIGHTTIME/DAYTIME INDICES
            nighttimeIndex = (data.time == 0) & (data.PPFD < 4);
            ntdata = data(nighttimeIndex, :);
            dtdata = data(~nighttimeIndex, :);
            
            %% STEP 10: CREATE DIAGNOSTIC FIGURES (for CRK)
            if strcmp(sitename, 'CRK')
                set(0, 'DefaultFigureVisible', 'off');
                
                % Convert -9999 to NaN for NEE columns
                data.NEE_orig(data.NEE_orig == -9999) = NaN;
                data.NEE_fall(data.NEE_fall == -9999) = NaN;
                
                % FIGURE 1: NT Respiration Model (5 panels)
                create_nt_respiration_figure(data, ntdata, sitename, yyear, config);
                
                % FIGURE 2: DT Respiration Model (5 panels)
                create_dt_respiration_figure(data, ntdata, sitename, yyear, config);
                
                % FIGURE 3: DT vs NT Comparison - Reco
                create_comparison_figure(data, sitename, yyear, config, 'Reco');
                
                % FIGURE 4: DT vs NT Comparison - GPP
                create_comparison_figure(data, sitename, yyear, config, 'GPP');
                
                set(0, 'DefaultFigureVisible', 'on');
            end
            
            %% STEP 11: PROCESS DAILY AND MONTHLY DATA
            smallertable = data(:, {'dt', 'GPP_f', 'Reco', 'NEE_f', 'VPD', 'PPFD', 'Reco_org', 'R_ref','Ustar','GPP_NT_f','Reco_NT_f'});
            data_tt = table2timetable(smallertable);
            
            % multiplier 
            12 * 10^(-6) *60*60*24*days
            % Convert flux units
            % 30 mins gC/m2/30 min 
            % Conversion multiplier for flux units
%multiplier = 12*60*60*24*10^(-6);
%secs_per_interval = 60 * 30; % 30-minute intervals
% 10*-6 converts micro mol to mol
% * 12.. 1 mole of CO2 has 1 mol C atoms 
% Atomic weigt of C = 12g/mol
% so 1 mol CO2 ocntains 12 g carbon 
% Below calculation gives gC/m2/30min
            data_tt.GPP_f = data_tt.GPP_f * 12 * 10^(-6) * secs_per_interval;
            data_tt.Reco = data_tt.Reco * 12 * 10^(-6) * secs_per_interval;
            data_tt.NEE_f = data_tt.NEE_f * 12 * 10^(-6) * secs_per_interval;
            data_tt.GPP_NT_f = data_tt.GPP_NT_f * 12 * 10^(-6) * secs_per_interval;
            data_tt.Reco_NT_f = data_tt.Reco_NT_f * 12 * 10^(-6) * secs_per_interval;
            
            % Daily aggregation
            % gC/m2/day
            daily_fluxes = retime(data_tt(:, {'GPP_f', 'Reco', 'NEE_f','GPP_NT_f','Reco_NT_f'}), 'daily', 'sum');
            daily_others = retime(data_tt(:, {'VPD', 'PPFD', 'Reco_org', 'R_ref','Ustar'}), 'daily', 'mean');
            daily_others.Reco_org = daily_others.Reco_org * 12 * 10^(-6) * secs_per_interval * 48;
            daily_others.R_ref = daily_others.R_ref * 12 * 10^(-6) * secs_per_interval * 48;
            daily = synchronize(daily_fluxes, daily_others);
            
            % Change zeros to NaN for flux variables
            daily.GPP_f(daily.GPP_f == 0) = NaN;
            daily.Reco(daily.Reco == 0) = NaN;
            daily.NEE_f(daily.NEE_f == 0) = NaN;
            
            %% STEP 12: MONTHLY AGGREGATION WITH GAP-FILLING
            daily_with_dt = daily;
            daily_with_dt.has_data = ~isnan(daily.NEE_f);
            month_days_available = retime(daily_with_dt(:, 'has_data'), 'monthly', 'sum');
            
            year_months = datetime(yyear, 1:12, 1);
            days_in_month = eomday(yyear, 1:12)';
            expected_days = array2table(days_in_month, 'VariableNames', {'days_in_month'});
            expected_days.dt = year_months';
            expected_days = table2timetable(expected_days);
            
            monthly_coverage = synchronize(month_days_available, expected_days);
            monthly_coverage.coverage = monthly_coverage.has_data ./ monthly_coverage.days_in_month;
            
            monthly_fluxes = retime(daily(:, {'GPP_f', 'Reco', 'NEE_f','GPP_NT_f','Reco_NT_f'}), 'monthly', 'sum');
            monthly_others = retime(daily(:, {'VPD', 'PPFD', 'Reco_org', 'R_ref'}), 'monthly', 'mean');
            monthly = synchronize(monthly_fluxes, monthly_others, monthly_coverage);
            
            monthly.Reco_org = monthly.Reco_org .* monthly.days_in_month;
            monthly.R_ref = monthly.R_ref .* monthly.days_in_month;
            
            % Gap-fill monthly sums
            monthly.GPP_f_filled = monthly.GPP_NT_f;
            monthly.Reco_filled = monthly.Reco_NT_f;
            monthly.NEE_f_filled = monthly.NEE_f;
            
            fill_idx = monthly.coverage >= 0.2 & monthly.coverage < 1;
            if any(fill_idx)
                monthly.GPP_f_filled(fill_idx) = monthly.GPP_NT_f(fill_idx) ./ monthly.coverage(fill_idx);
                monthly.Reco_filled(fill_idx) = monthly.Reco_NT_f(fill_idx) ./ monthly.coverage(fill_idx);
                monthly.NEE_f_filled(fill_idx) = monthly.NEE_f(fill_idx) ./ monthly.coverage(fill_idx);
            end
            
            monthly.Month = month(monthly.dt);
            monthly.Year = year(monthly.dt);
            
            %% STEP 13: SPECIAL CASE - NC2 2012 MANUAL VALUES
            if yyear == 2012 && strcmp(sitename, 'NC2')
                monthly = apply_nc2_2012_corrections(monthly, months);
            end
            
            %% STEP 13B: SPECIAL CASE - CRK 2023 MANUAL VALUES
            if yyear == 2023 && strcmp(sitename, 'CRK')
                monthly = apply_crk_2023_corrections(monthly, months);
            end
            
            %% STEP 14: LOAD AND PROCESS METEOROLOGICAL DATA
            directory = [config.input_dir, '\Annual data\', sitename];
            filename = sprintf('%s_%d.csv', sitename, yyear);
            fullFilePath = fullfile(directory, filename);
            
            if ~isfile(fullFilePath)
                fprintf('⚠ (met data not found) ');
                met_table_daily = [];
                met_table_monthly = [];
            else
                all_data = readtable(fullFilePath);
                
                all_data.dt = datetime(all_data.Year, 1, 0) + caldays(all_data.doy) + hours(all_data.hours);
                
                % Replace -9999 with NaN
                varNames = all_data.Properties.VariableNames;
                for i = 1:numel(varNames)
                    variable = all_data.(varNames{i});
                    if isnumeric(variable)
                        variable(variable == -9999) = NaN;
                    end
                    all_data.(varNames{i}) = variable;
                end
                
                % Site-specific met processing
                switch sitename
                    case 'NC2'
                        [met_table_daily, met_table_monthly] = process_nc2_met(all_data, data);
                    case 'CRK'
                        [met_table_daily, met_table_monthly] = process_crk_met(all_data, data, yyear);
                    case 'CST'
                        [met_table_daily, met_table_monthly] = process_cst_met(all_data, data);
                    case 'GA'
                        [met_table_daily, met_table_monthly] = process_ga_met(all_data, data);
                end
            end
            
            %% STEP 15: PROCESS SOIL RESPIRATION DATA
            switch sitename
                case 'NC2'
                    merged_data_monthly = process_nc2_sr(monthly, met_table_monthly, yyear, config);
                case 'GA'
                    merged_data_monthly = process_ga_sr(monthly, met_table_monthly, config);
                case 'CRK'
                    merged_data_monthly = process_crk_sr(monthly, met_table_monthly, yyear, config);
                case 'CST'
                    if ~isempty(met_table_monthly)
                        merged_data_monthly = synchronize(monthly, met_table_monthly, 'union');
                    else
                        merged_data_monthly = monthly;
                    end
            end
            
            %% STEP 16: SAVE OUTPUT FILES
            directory = [config.output_dir, '\TIMBCA\Monthly\Seasonal\', sitename];
            if ~exist(directory, 'dir')
                mkdir(directory);
            end
            
            % Save monthly data
            filename = sprintf('%s_%d_monthly.csv', sitename, yyear);
            fullFilePath = fullfile(directory, filename);
            writetimetable(merged_data_monthly, fullFilePath);
            
            % Save daily data
            if ~isempty(met_table_daily)
                merged_data_daily = synchronize(daily, met_table_daily, 'union');
            else
                merged_data_daily = daily;
            end
            filename = sprintf('%s_%d_daily.csv', sitename, yyear);
            fullFilePath = fullfile(directory, filename);
            writetimetable(merged_data_daily, fullFilePath);
            
            fprintf('✓ COMPLETE\n');
            total_processed = total_processed + 1;
            
        catch ME
            total_errors = total_errors + 1;
            continue;
        end
    end
end

%% FINAL SUMMARY
fprintf('\n╔═══════════════════════════════════════════════════════════════╗\n');
fprintf('║  BATCH PROCESSING COMPLETE                                    ║\n');
fprintf('╚═══════════════════════════════════════════════════════════════╝\n\n');
fprintf('  Successfully processed:  %d site-year(s)\n', total_processed);
fprintf('  Skipped (not found):     %d site-year(s)\n', total_skipped);
fprintf('  Errors:                  %d site-year(s)\n\n', total_errors);

%% ═══════════════════════════════════════════════════════════════════════
%  HELPER FUNCTIONS
%  ═══════════════════════════════════════════════════════════════════════

function create_nt_respiration_figure(data, ntdata, sitename, yyear, config)
    % Create NT respiration model figure with 5 panels
    fig = figure('Position', [100, 100, 1800, 600], 'Color', 'w');
    colors = struct('measured',[0.2,0.6,0.2], 'modeled',[0.8,0.2,0.2], 'gpp',[0.6,0.2,0.8], 'fit',[0.2,0.2,0.2]);
    
    % Filter quality data
    if strcmp(sitename,'CRK') && yyear == 2024
        quality_data = ntdata(ntdata.QC_flag == 0 | ntdata.QC_flag == 1, :);
    else
        quality_data = ntdata(ntdata.QC_flag == 0, :);
    end
    valid_idx = ~isnan(quality_data.Reco_org) & ~isnan(quality_data.Reco);
    valid_data = quality_data(valid_idx, :);
    
    % Panel 1: Respiration time series
    subplot(1,5,1);
    scatter(data.dt, data.Reco_org, 30, 'o', 'MarkerFaceColor', colors.measured, 'MarkerEdgeColor', colors.measured, 'MarkerFaceAlpha',0.6);
    hold on;
    scatter(data.dt, data.Reco, 20, 'x', 'MarkerEdgeColor', colors.modeled, 'LineWidth',1.5);
    xlabel('Date','FontWeight','bold');
    ylabel('Respiration (μmol CO_2 m^{-2} s^{-1})');
    title('Respiration Time Series');
    legend({'Measured Reco','Modeled Reco'},'Location','best');
    datetick('x','mmm-yy','keepticks'); xtickangle(45);
    grid on; box on;
    
    % Panel 2: GPP
    subplot(1,5,2);
    scatter(data.dt, data.GPP_f, 25, 's', 'MarkerFaceColor', colors.gpp, 'MarkerEdgeColor', colors.gpp, 'MarkerFaceAlpha',0.6);
    xlabel('Date','FontWeight','bold');
    ylabel('GPP (μmol CO_2 m^{-2} s^{-1})');
    title('Gross Primary Production');
    legend('GPP','Location','best');
    datetick('x','mmm-yy','keepticks'); xtickangle(45);
    grid on; box on;
    
    % Panel 3: Model performance
    subplot(1,5,3);
    if ~isempty(valid_data)
        scatter(valid_data.Reco_org, valid_data.Reco, 40, 'o', 'MarkerFaceColor', colors.modeled, 'MarkerEdgeColor', colors.modeled, 'MarkerFaceAlpha',0.5);
        hold on;
        lm = fitlm(valid_data.Reco_org, valid_data.Reco);
        x_range = linspace(min(valid_data.Reco_org), max(valid_data.Reco_org), 100);
        plot(x_range, predict(lm, x_range'), '-', 'Color', colors.fit, 'LineWidth', 2);
        axis_lims = [min([valid_data.Reco_org; valid_data.Reco]) max([valid_data.Reco_org; valid_data.Reco])];
        plot(axis_lims, axis_lims, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
        rmse = sqrt(mean((valid_data.Reco - valid_data.Reco_org).^2));
        text(0.05, 0.95, sprintf('R^2 = %.3f\nRMSE = %.3f', lm.Rsquared.Ordinary, rmse), 'Units','normalized', 'VerticalAlignment','top', 'BackgroundColor','w', 'EdgeColor','k');
        xlabel('Measured Reco (QC = 0)','FontWeight','bold');
        ylabel('Modeled Reco','FontWeight','bold');
        title('Model Performance');
        axis square; grid on; box on;
    end
    
    % Panel 4: NEE time series
    subplot(1,5,4);
    scatter(data.dt, data.NEE_orig, 20, 'o', 'MarkerFaceColor', colors.measured, 'MarkerEdgeColor', colors.measured, 'MarkerFaceAlpha', 0.6);
    hold on;
    scatter(data.dt, data.NEE_fall, 20, 'x', 'MarkerEdgeColor', colors.modeled, 'LineWidth', 1.5);
    legend("NEE (orig)", "NEE (modeled)", 'Location', 'best');
    xlabel('Date'); ylabel('NEE (μmol CO_2 m^{-2} s^{-1})');
    title('NEE Time Series');
    datetick('x','mmm-yy','keepticks'); xtickangle(45);
    grid on; box on;
    
    % Panel 5: NEE model performance
    subplot(1,5,5);
    scatter(data.NEE_orig, data.NEE_fall, 40, 'o', 'MarkerFaceColor', colors.modeled, 'MarkerEdgeColor', colors.modeled, 'MarkerFaceAlpha', 0.5);
    hold on;
    lm_nee = fitlm(data, 'NEE_fall~NEE_orig');
    plot(lm_nee, 'color', colors.fit);
    xlabel('Measured NEE'); ylabel('Modeled NEE');
    title('NEE Model Performance');
    xLimits = xlim; yLimits = ylim;
    text(xLimits(1) + 0.05*range(xLimits), yLimits(2) - 0.2*range(yLimits), sprintf('R^2 = %.3f', lm_nee.Rsquared.ordinary), 'BackgroundColor','w','EdgeColor','k');
    grid on; box on;
    
    sgtitle(sprintf('NT Respiration Model - %s %d', sitename, yyear), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save
    filename = sprintf('NT_Respiration_NEE_Combined_%s_%d.png', sitename, yyear);
    saveas(fig, fullfile(config.fig_dir, filename));
    close(fig);
end

function create_dt_respiration_figure(data, ntdata, sitename, yyear, config)
    % Create DT respiration model figure with 5 panels
    fig = figure('Position', [100, 100, 1800, 600], 'Color', 'w');
    colors = struct('measured',[0.2,0.6,0.2], 'modeled',[0.8,0.2,0.2], 'gpp',[0.6,0.2,0.8], 'fit',[0.2,0.2,0.2]);
    
    quality_data = ntdata(ntdata.QC_flag == 0, :);
    valid_idx = ~isnan(quality_data.Reco_org) & ~isnan(quality_data.Reco_DT);
    valid_data = quality_data(valid_idx, :);
    
    % Panel 1: Respiration time series
    subplot(1,5,1);
    scatter(data.dt, data.Reco_org, 30, 'o', 'MarkerFaceColor', colors.measured, 'MarkerEdgeColor', colors.measured, 'MarkerFaceAlpha',0.6);
    hold on;
    scatter(data.dt, data.Reco_DT, 20, 'x', 'MarkerEdgeColor', colors.modeled, 'LineWidth',1.5);
    xlabel('Date','FontWeight','bold');
    ylabel('Respiration (μmol CO_2 m^{-2} s^{-1})');
    title('Respiration Time Series');
    legend({'Measured Reco','Modeled Reco'},'Location','best');
    datetick('x','mmm-yy','keepticks'); xtickangle(45);
    grid on; box on;
    
    % Panel 2: GPP
    subplot(1,5,2);
    scatter(data.dt, data.GPP_DT, 25, 's', 'MarkerFaceColor', colors.gpp, 'MarkerEdgeColor', colors.gpp, 'MarkerFaceAlpha',0.6);
    xlabel('Date','FontWeight','bold');
    ylabel('GPP (μmol CO_2 m^{-2} s^{-1})');
    title('Gross Primary Production');
    legend('GPP','Location','best');
    datetick('x','mmm-yy','keepticks'); xtickangle(45);
    grid on; box on;
    
    % Panel 3: Model performance
    subplot(1,5,3);
    if ~isempty(valid_data)
        scatter(valid_data.Reco_org, valid_data.Reco_DT, 40, 'o', 'MarkerFaceColor', colors.modeled, 'MarkerEdgeColor', colors.modeled, 'MarkerFaceAlpha',0.5);
        hold on;
        lm = fitlm(valid_data.Reco_org, valid_data.Reco_DT);
        x_range = linspace(min(valid_data.Reco_org), max(valid_data.Reco_org), 100);
        plot(x_range, predict(lm, x_range'), '-', 'Color', colors.fit, 'LineWidth', 2);
        axis_lims = [min([valid_data.Reco_org; valid_data.Reco_DT]) max([valid_data.Reco_org; valid_data.Reco_DT])];
        plot(axis_lims, axis_lims, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
        rmse = sqrt(mean((valid_data.Reco_DT - valid_data.Reco_org).^2));
        text(0.05, 0.95, sprintf('R^2 = %.3f\nRMSE = %.3f', lm.Rsquared.Ordinary, rmse), 'Units','normalized', 'VerticalAlignment','top', 'BackgroundColor','w', 'EdgeColor','k');
        xlabel('Measured Reco (QC = 0)','FontWeight','bold');
        ylabel('Modeled Reco','FontWeight','bold');
        title('Model Performance');
        axis square; grid on; box on;
    end
    
    % Panel 4: NEE time series
    subplot(1,5,4);
    scatter(data.dt, data.NEE_orig, 20, 'o', 'MarkerFaceColor', colors.measured, 'MarkerEdgeColor', colors.measured, 'MarkerFaceAlpha', 0.6);
    hold on;
    scatter(data.dt, data.NEE_fall, 20, 'x', 'MarkerEdgeColor', colors.modeled, 'LineWidth', 1.5);
    legend("NEE (orig)", "NEE (modeled)", 'Location', 'best');
    xlabel('Date'); ylabel('NEE (μmol CO_2 m^{-2} s^{-1})');
    title('NEE Time Series');
    datetick('x','mmm-yy','keepticks'); xtickangle(45);
    grid on; box on;
    
    % Panel 5: NEE model performance
    subplot(1,5,5);
    scatter(data.NEE_orig, data.NEE_fall, 40, 'o', 'MarkerFaceColor', colors.modeled, 'MarkerEdgeColor', colors.modeled, 'MarkerFaceAlpha', 0.5);
    hold on;
    lm_nee = fitlm(data, 'NEE_fall~NEE_orig');
    plot(lm_nee, 'color', colors.fit);
    xlabel('Measured NEE'); ylabel('Modeled NEE');
    title('NEE Model Performance');
    xLimits = xlim; yLimits = ylim;
    text(xLimits(1) + 0.05*range(xLimits), yLimits(2) - 0.2*range(yLimits), sprintf('R^2 = %.3f', lm_nee.Rsquared.ordinary), 'BackgroundColor','w','EdgeColor','k');
    grid on; box on;
    
    sgtitle(sprintf('DT Respiration Model - %s %d', sitename, yyear), 'FontSize', 16, 'FontWeight', 'bold');
    
    % Save
    filename = sprintf('DT_Respiration_NEE_Combined_%s_%d.png', sitename, yyear);
    saveas(fig, fullfile(config.fig_dir, filename));
    close(fig);
end

function create_comparison_figure(data, sitename, yyear, config, flux_type)
    % Create DT vs NT comparison figure
    fig = figure('Position', [100, 100, 1200, 800]);
    
    if strcmp(flux_type, 'Reco')
        scatter(data.dt, data.Reco_DT, 20, 'o', 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0.3010, 0.7450, 0.9330], 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Reco DT');
        hold on;
        scatter(data.dt, data.Reco, 20, 'o', 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor', [0.6780, 0.9210, 0.4430], 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Reco NT');
        ylabel('Reco (μmol CO_2 m^{-2} s^{-1})');
        title(sprintf('Daytime vs Nighttime RECO - %s %d', sitename, yyear));
        filename = sprintf('daytime_vs_nighttime_RECO_%s_%d.png', sitename, yyear);
    else % GPP
        scatter(data.dt, data.GPP_DT, 20, 'o', 'MarkerEdgeColor', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0.3010, 0.7450, 0.9330], 'MarkerFaceAlpha', 0.3, 'DisplayName', 'GPP DT');
        hold on;
        scatter(data.dt, data.GPP_f, 20, 'o', 'MarkerEdgeColor', [0.4660, 0.6740, 0.1880], 'MarkerFaceColor', [0.6780, 0.9210, 0.4430], 'MarkerFaceAlpha', 0.3, 'DisplayName', 'GPP NT');
        ylabel('GPP (μmol CO_2 m^{-2} s^{-1})');
        title(sprintf('Daytime vs Nighttime GPP - %s %d', sitename, yyear));
        filename = sprintf('daytime_vs_nighttime_GPP_%s_%d.png', sitename, yyear);
    end
    
    xlabel('Date');
    legend('Location', 'best');
    grid on;
    
    saveas(fig, fullfile(config.fig_dir, filename));
    close(fig);
end

function monthly = apply_nc2_2012_corrections(monthly, months)
    % Apply manual corrections for NC2 2012
    months_to_update = [1, 2, 3, 4];
    manual_values = struct();
    manual_values(1).GPP = 99.42; manual_values(1).Reco = 67; manual_values(1).NEE = -32;
    manual_values(2).GPP = 105; manual_values(2).Reco = 60; manual_values(2).NEE = -45;
    manual_values(3).GPP = 166; manual_values(3).Reco = 79; manual_values(3).NEE = -87;
    manual_values(4).GPP = 252; manual_values(4).Reco = 176; manual_values(4).NEE = -77;
    
    monthly_date_only = dateshift(monthly.dt, 'start', 'day');
    
    for month_idx = months_to_update
        specific_date = datetime(['01-', char(months(month_idx)), '-2012'], 'Format', 'dd-MMM-yyyy');
        if any(monthly_date_only == specific_date)
            monthly.GPP_f(monthly_date_only == specific_date) = manual_values(month_idx).GPP;
            monthly.Reco(monthly_date_only == specific_date) = manual_values(month_idx).Reco;
            monthly.NEE_f(monthly_date_only == specific_date) = manual_values(month_idx).NEE;
            monthly.GPP_NT_f(monthly_date_only == specific_date) = manual_values(month_idx).GPP;
            monthly.Reco_NT_f(monthly_date_only == specific_date) = manual_values(month_idx).Reco;
        end
    end
end

function monthly = apply_crk_2023_corrections(monthly, months)
    % Apply manual corrections for CRK 2023 (Sept, Oct, Nov, Dec)
    months_to_update = [9, 10, 11, 12];  % Sept through Dec
    
    % Define manual values for each month (from your image)
    manual_values = struct();
    
    % September
    manual_values(9).GPP_f = NaN;
    manual_values(9).Reco = NaN;
    manual_values(9).NEE_f = 14.67722;
    manual_values(9).GPP_NT_f = 126.816;
    manual_values(9).Reco_NT_f = 141.4932;
    
    % October
    manual_values(10).GPP_f = NaN;
    manual_values(10).Reco = NaN;
    manual_values(10).NEE_f = -4.83;
    manual_values(10).GPP_NT_f = 125.41;
    manual_values(10).Reco_NT_f = 120.58;
    
    % November
    manual_values(11).GPP_f = NaN;
    manual_values(11).Reco = NaN;
    manual_values(11).NEE_f = 46.01;
    manual_values(11).GPP_NT_f = 103.94;
    manual_values(11).Reco_NT_f = 149.95;
    
    % December
    manual_values(12).GPP_f = NaN;
    manual_values(12).Reco = NaN;
    manual_values(12).NEE_f = 39.744;
    manual_values(12).GPP_NT_f = 75.796;
    manual_values(12).Reco_NT_f = 115.55;
    
    monthly_date_only = dateshift(monthly.dt, 'start', 'day');
    
    for month_idx = months_to_update
        specific_date = datetime(['01-', char(months(month_idx)), '-2023'], 'Format', 'dd-MMM-yyyy');
        if any(monthly_date_only == specific_date)
            monthly.GPP_f(monthly_date_only == specific_date) = manual_values(month_idx).GPP_f;
            monthly.Reco(monthly_date_only == specific_date) = manual_values(month_idx).Reco;
            monthly.NEE_f(monthly_date_only == specific_date) = manual_values(month_idx).NEE_f;
            monthly.GPP_NT_f(monthly_date_only == specific_date) = manual_values(month_idx).GPP_NT_f;
            monthly.Reco_NT_f(monthly_date_only == specific_date) = manual_values(month_idx).Reco_NT_f;
        end
    end
end

function [met_table_daily, met_table_monthly] = process_nc2_met(all_data, data)
    % Extract met data columns
    met_table = all_data(:, {'dt','SWC_1_1_1','P_1_1_1','VPD_PI_1_1_1', 'TS_1_2_1'});
    
    % Create a table from data with Tair_f and Tsoil
    temp_data = table(data.dt, data.Tair_f, data.Tsoil, 'VariableNames', {'dt', 'Tair_f', 'Tsoil'});
    
    % Merge based on datetime to handle dimension mismatches
    met_table = outerjoin(met_table, temp_data, 'Keys', 'dt', 'MergeKeys', true);
    
    % Convert to timetable and aggregate
    met_table_tt = table2timetable(met_table);
    met_table_daily = retime(met_table_tt, 'daily', 'mean');
    met_table_monthly = retime(met_table_daily, 'monthly', 'mean');
end

function [met_table_daily, met_table_monthly] = process_crk_met(all_data, data, year)
    met_table = all_data(:, {'dt', 'SWC_1_1_1', 'SWC_1_2_1', 'P_1_1_1', 'VPD_1_1_1', ...
                             'TS_1_2_1', 'TS_1_1_1', 'TA_1_1_1', 'NETRAD_1_1_1', 'PPFD_IN_1_1_1'});
    
    % Clean data
    varNames = met_table.Properties.VariableNames;
    for i = 1:numel(varNames)
        variable = met_table.(varNames{i});
        if isnumeric(variable)
            variable(variable == -9999) = NaN;
        end
        met_table.(varNames{i}) = variable;
    end
    
    % Handle missing row
    met_table(end+1, :) = met_table(end-1, :);
    
    % Set date for last row based on year
    if year == 2022
        met_table.dt(end) = datetime('1-Jan-2023 00:00:00');
    elseif year == 2023
        met_table.dt(end) = datetime('1-Jan-2024 00:00:00');
    elseif year == 2024
        met_table.dt(end) = datetime('1-Jan-2025 00:00:00');
    end
    
    % Update values for last row
    met_table.SWC_1_1_1(end) = 22.4000;
    met_table.P_1_1_1(end) = 0;
    
    % Add data from other source
    met_table.SWC_5cm_f = data.SWC_5cm_f;
    met_table.SWC_20cm_f = data.SWC_20cm_f;
    
    % Add temperature data
    met_table.Tair_f = data.Tair_f;
    met_table.Tsoil = data.Tsoil;
    
    % Process to daily and monthly
    met_table_tt = table2timetable(met_table);
    
    % Mean values for most meteorological variables
    mean_tt = retime(met_table_tt(:, {'SWC_1_1_1', 'SWC_1_2_1', 'Tair_f', 'TS_1_1_1', 'TS_1_2_1', ...
                                      'Tsoil', 'VPD_1_1_1', 'TA_1_1_1', 'SWC_5cm_f', 'SWC_20cm_f', ...
                                      'NETRAD_1_1_1', 'PPFD_IN_1_1_1'}), 'daily', 'mean');
    mean_tt_monthly = retime(mean_tt, 'monthly', 'mean');
    
    % Sum values for precipitation
    p_tt = retime(met_table_tt(:, {'P_1_1_1'}), 'daily', 'sum');
    p_tt_monthly = retime(p_tt, 'monthly', 'sum');
    
    % Merge tables
    met_table_daily = synchronize(mean_tt, p_tt, 'intersection');
    met_table_monthly = synchronize(mean_tt_monthly, p_tt_monthly, 'intersection');
end

function [met_table_daily, met_table_monthly] = process_cst_met(all_data, data)
    % Extract met data columns
    met_table = all_data(:, {'dt','SWC_1_1_1','P'});
    
    % Create a table from data with Tair and Tsoil
    temp_data = table(data.dt, data.Tair, data.Tsoil, 'VariableNames', {'dt', 'Tair_f', 'Tsoil'});
    
    % Merge based on datetime to handle dimension mismatches
    met_table = outerjoin(met_table, temp_data, 'Keys', 'dt', 'MergeKeys', true);
    
    % Convert to timetable and aggregate
    met_table_tt = table2timetable(met_table);
    met_table_daily = retime(met_table_tt, 'daily', 'mean');
    met_table_monthly = retime(met_table_daily, 'monthly', 'mean');
end

function [met_table_daily, met_table_monthly] = process_ga_met(all_data, data)
    % Extract met data columns
    met_table = all_data(:, {'dt','SWC','P_RAIN'});
    
    % Create a table from data with Tair and Tsoil
    temp_data = table(data.dt, data.Tair, data.Tsoil, 'VariableNames', {'dt', 'Tair_f', 'Tsoil'});
    
    % Merge based on datetime to handle dimension mismatches
    met_table = outerjoin(met_table, temp_data, 'Keys', 'dt', 'MergeKeys', true);
    
    % Convert to timetable and aggregate
    met_table_tt = table2timetable(met_table);
    met_table_daily = retime(met_table_tt, 'daily', 'mean');
    met_table_monthly = retime(met_table_daily, 'monthly', 'mean');
end

function merged_data_monthly = process_nc2_sr(monthly, met_table_monthly, yyear, config)
    sr_file = [config.output_dir, '\TIMBCA\Monthly\SR\NC2\NC2_monthly_SR.csv'];
    if isfile(sr_file)
        sr = readtable(sr_file);
        index = (sr.Year == yyear);
        sr = sr(index, :);
        if ~isempty(sr)
            dateStrings = strcat(string(sr.Year), "-", sr.Month);
            sr.dt = datetime(dateStrings, 'InputFormat', 'yyyy-MMM');
            monthly_sr = table2timetable(sr);
            merged_data_monthly = synchronize(monthly, met_table_monthly, monthly_sr, 'union');
        else
            merged_data_monthly = synchronize(monthly, met_table_monthly, 'union');
        end
    else
        merged_data_monthly = synchronize(monthly, met_table_monthly, 'union');
    end
end

function merged_data_monthly = process_crk_sr(monthly, met_table_monthly, yyear, config)
    sr_file = [config.output_dir, '\TIMBCA\Monthly\SR\CRK\CRK_monthly_SR.csv'];
    if isfile(sr_file)
        sr = readtable(sr_file);
        index = (sr.Year == yyear);
        sr = sr(index, :);
        if ~isempty(sr)
            dateStrings = strcat(string(sr.Year), "-", sr.Month);
            sr.dt = datetime(dateStrings, 'InputFormat', 'yyyy-MMM');
            monthly_sr = table2timetable(sr);
            merged_data_monthly = synchronize(monthly, met_table_monthly, monthly_sr, 'union');
        else
            merged_data_monthly = synchronize(monthly, met_table_monthly, 'union');
        end
    else
        merged_data_monthly = synchronize(monthly, met_table_monthly, 'union');
    end
end

function merged_data_monthly = process_ga_sr(monthly, met_table_monthly, config)
    merged_data_monthly = synchronize(monthly, met_table_monthly, 'union');
end