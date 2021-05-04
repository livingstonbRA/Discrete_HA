%% ONE ASSET DISCRETE TIME HA MODEL
% This is the main script for this code repository.
% See the readme for details.

clear;
close all;

%% ------------------------------------------------------------------------
% SET OPTIONS
% -------------------------------------------------------------------------
% options
runopts.calibrate = true; % wrap code in nonlinear solver
runopts.fast = false; % very small asset and income grids for testing
runopts.Simulate = false; % also solve distribution via simulation
runopts.MakePlots = true;
runopts.MPCs = true;
runopts.MPCs_news = false;
runopts.MPCs_loan_and_loss = false;
runopts.DeterministicMPCs = true; % must be on if decompositions are needed
runopts.SaveOutput = true;

% name of parameters script in code/+params directory
runopts.mode = 'parameters'; % 'parameters'

% select experiment (ignored when run on server)
runopts.name_to_run = ''; % ''
runopts.number = [2]; % []

%% ------------------------------------------------------------------------
% HOUSEKEEPING, DO NOT CHANGE
% -------------------------------------------------------------------------
% Get task id if running on server
server_array_id = str2num(getenv('SLURM_ARRAY_TASK_ID'));
running_on_server = ~isempty(server_array_id);
if running_on_server
    runopts.number = server_array_id;
end

% Path for .mat output file
matname = sprintf('variables%d.mat', runopts.number);
runopts.savematpath = fullfile('output', matname);
xlxname = sprintf('table%d.xlsx', runopts.number);
runopts.savexlxpath = fullfile('output', xlxname);

% Directories
warning('off', 'MATLAB:MKDIR:DirectoryExists')
mkdir('output')
mkdir('temp')

if exist(runopts.savematpath, 'file') == 2
    % Delete old results
    delete runopts.savematpath;
end

addpath('code');
addpath(fullfile('code', 'aux_lib'));

%% ------------------------------------------------------------------------
% LOAD PARAMETERS
% -------------------------------------------------------------------------
params = setup.(runopts.mode)(runopts);

%% ------------------------------------------------------------------------
% CALL MAIN FUNCTION
% -------------------------------------------------------------------------
fprintf('\nParameterization "%s" was chosen.\n', params.name)

if params.calibrating
    % Calibrate with nonlinear solver
    disp('Beginning model calibration...')
    calibrator = params.calibrator;
    options = optimoptions(@lsqnonlin, 'MaxIterations', params.calibrate_maxiter,...
            'FunctionTolerance', params.calibrate_tol);
    solver_args = params.calibrator.get_args();
    calibrated_params = lsqnonlin(params.calibrator.solver_handle,...
        solver_args{:}, options);

    params.calibrating = false;
end

% Final run
results = main(params);
fprintf('Finished parameterization %s\n', params.name)

%% ------------------------------------------------------------------------
% CREATE TABLE OF RESULTS
% -------------------------------------------------------------------------
table_out = tables.OutputTable(params, results.stats)
writetable(table_out, runopts.savexlxpath, 'WriteRowNames', true);

%% ------------------------------------------------------------------------
% FIGURES
% -------------------------------------------------------------------------
if runopts.MakePlots
    mkdir('output/figures')
    if strcmp(params.name, 'Quarterly')
        % Quarterline baseline
        statistics.DHAPlots.consumption_wealth_overlay(results.stats, params);
        figpath = fullfile('output', 'baseline_quarterly_c_x_overlay.jpg');
        saveas(gcf, figpath)
        
        statistics.DHAPlots.mpc_function(results.stats, 'mpc_type', 'period');
        figpath = fullfile('output', 'baseline_quarterly_period_mpcs.jpg');
        saveas(gcf, figpath)
        
        statistics.DHAPlots.mpc_function(results.stats, 'mpc_type', 'cumulative');
        figpath = fullfile('output', 'baseline_quarterly_cumulative_mpcs.jpg');
        saveas(gcf, figpath)
    elseif params.nbeta > 1
        statistics.DHAPlots.consumption_wealth_overlay(results.stats, params,...
            'curve_variable', 'beta');
        figpath = fullfile('output', 'fixed_beta_heterogeneity_c_x_overlay.jpg');
        saveas(gcf, figpath)
    end
end
        