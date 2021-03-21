function [params, all_names] = parameters(runopts)
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    import solver.DHACalibrator
    
    scf = setup.scf2019struct();

    dollars = [-1, -500, -5000, 1, 500, 5000];
    shared_params.annual_inc_dollars = scf.quarterly_earnings * 4;
    shared_params.shocks = dollars ./ shared_params.annual_inc_dollars;

    shared_params.xgrid_par = 0.1;
    shared_params.xgrid_term1wt = 0.01;
    shared_params.xgrid_term1curv = 0.5;
    shared_params.xmax = 500;

    shared_params.shocks_labels = {};
    for ishock = 1:6
        val = dollars(ishock);
        if val < 0
            shared_params.shocks_labels{ishock} = sprintf('-$%g', abs(val));
        else
            shared_params.shocks_labels{ishock} = sprintf('$%g', abs(val));
        end
    end

   income_path = 'input/income_quarterly_b.mat';

   % A parameterization
   params(1) = setup.Params(4, 'Quarterly Model', income_path);

   %----------------------------------------------------------------------
    % CALL METHODS/CHANGE SELECTED PARAMETERS, DO NOT CHANGE
    %----------------------------------------------------------------------
    params.set_index();

    % get list of all names
    all_names = cell2table({params.name}');
    
    % select by number if there is one, otherwise select by names,
    % otherwise use all
    if numel(runopts.number) == 1
        params = params(runopts.number);
    elseif numel(runopts.number) > 1
        error('runopts.number must have 1 or zero elements')
    else
        params = setup.Params.select_by_names(params, runopts.name_to_run);
    end

    params.set_run_parameters(runopts);
    params.make_adjustments();

    %----------------------------------------------------------------------
    % ATTACH CALIBRATOR
    %----------------------------------------------------------------------
    if params.calibrate
        calibration = calibrations(params.index);
        calibrator = DHACalibrator(params, calibration.variables,...
            calibration.target_names, calibration.target_values);
        heterogeneity = setup.Prefs_R_Heterogeneity(params);

        if (params.nbeta > 1) && isequal(heterogeneity.ztrans, eye(params.nbeta))
            new_betaH = params.betaH - max(heterogeneity.betagrid0);
            params.set("betaH", new_betaH, true);
        end

        beta_bounds = [params.betaL, params.betaH];
        calibrator.set_param_bounds(beta_bounds);
        calibrator.set_handle(params);
        params.set("calibrator", calibrator, true);
    end
end