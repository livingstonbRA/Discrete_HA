function [params, all_names] = parameters(runopts)
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    import solver.DHACalibrator
    
    %% SET SHARED PARAMETERS
    shared_params = struct();
    shared_params.r = 0.01;
    shared_params.Bequests = false;

    % Statistics from the 2019 SCF
    scf = setup.scf2019struct();
    shared_params.numeraire_in_dollars = scf.annual_earnings;

    % Income processes
    quarterly_b_params = shared_params;
    quarterly_b_params.IncomeProcess = 'input/income_quarterly_b_contyT.mat';
    quarterly_b_params.sd_logyT = sqrt(0.6376);
    quarterly_b_params.lambdaT = 0.25;

    quarterly_c_params = shared_params;
    quarterly_c_params.IncomeProcess = 'input/income_quarterly_c_contyT.mat';
    quarterly_c_params.sd_logyT = sqrt(1.6243);
    quarterly_c_params.lambdaT = 0.0727;

    quarterly_a_params = shared_params;
    quarterly_a_params.sd_logyT = sqrt(0.2087);
    quarterly_a_params.sd_logyP = sqrt(0.01080);
    quarterly_a_params.rho_logyP = 0.9881;
    quarterly_a_params.lambdaT = 1;

    annual_params = shared_params;
    annual_params.sd_logyT = sqrt(0.0494);
    annual_params.sd_logyP = sqrt(0.0422);
    annual_params.rho_logyP = 0.9525;
    annual_params.lambdaT = 1;
    
    %% SET CALIBRATIONS
    % Ignored if not calibrating
    calibrations = struct();
    calibrations.variables = {'beta0'};
    % calibrations.target_names = {'median_a'};
    % calibrations.target_values = [scf.median_totw];
    calibrations.target_names = {'mean_a'};
    calibrations.target_values = [4.1];

    for ii = 2:999
        calibrations(ii) = calibrations(1);
    end

    %% PARAMETERIZATIONS

    %----------------------------------------------------------------------
    % BASELINES
    %----------------------------------------------------------------------
    
    % Annual
    params(1) = setup.Params(1, 'Annual', annual_params);
    params(1).beta0 = 0.984108034755346;
    params(1).group = {'Baseline', 'A1'};
    params(1).descr = 'Baseline';
     
    % Quarterly
    params(end+1) = setup.Params(4, 'Quarterly', quarterly_b_params);
    params(end).beta0 = 0.984363510593659;
    params(end).group = {'Baseline', 'Q1a', 'Q1b', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6'};
    params(end).descr = 'Baseline';
    params(end).tex_header = 'Baseline';

    %----------------------------------------------------------------------
    % PART 2, DIFFERENT ASSUMPTIONS
    %----------------------------------------------------------------------
    ifreq = 4;

    % Liquid wealth calibration, mean assets = 2.25
    name = sprintf('E[a] = %g', 2.25);
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).group = {'Q1a'};
    params(end).descr = 'Calibration to liquid wealth, E[a] = 2.25';
    params(end).tex_header = 'E[a]';
    params(end).tex_header_values = {struct('value', 2.25)};
    n = numel(params);
    calibrations(n).target_names = {'mean_a'};
    calibrations(n).target_values = [2.25];
    
    % Total wealth calibration, mean assets = 9.4
    name = sprintf('E[a] = %g', 9.4);
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).group = {'Q1a'};
    params(end).descr = {'Calibration to total wealth, E[a] = 9.4'};
    params(end).tex_header = 'E[a]';
    params(end).tex_header_values = {struct('value', 9.4)};
    n = numel(params);
    calibrations(n).target_names = {'mean_a'};
    calibrations(n).target_values = [9.4];
    
    % Liquid wealth calibration, median assets = 0.05, 0.5, 1.0
    for mw = [0.05, 0.5, 1.0]
        name = sprintf('median(a) = %g', mw);
        params(end+1) = setup.Params(4, name, quarterly_b_params);
        params(end).group = {'Q1a'};
        params(end).descr = {sprintf('Calibration to liquid wealth, median(a) = %g', mw)};
        params(end).tex_header = 'Median(a)';
        params(end).tex_header_values = {struct('value', mw)};
        n = numel(params);
        calibrations(n).target_names = {'median_a'};
        calibrations(n).target_values = [mw];
    end

    % no death
    name = 'No Death';
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).group = {'Q1b'};
    params(end).dieprob = 0;
    params(end).beta0 = 0.975363510593659;
    params(end).tex_header = 'No Death';

    % % no bequests
    % name = 'No Bequests';
    % params(end+1) = setup.Params(4, name, quarterly_b_params);
    % params(end).group = {'Q1b'};
    % params(end).Bequests = 0;
    % params(end).tex_header = 'No Bequests';

    % With bequests
    name = 'With Bequests';
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).group = {'Q1b'};
    params(end).Bequests = true;
    params(end).tex_header = 'With Bequests';

    % perfect annuities
    name = 'Annuities';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).group = {'Q1b'};
    params(end).annuities = true;
    params(end).betaH0 = - 5e-3;
    params(end).tex_header = 'Annuities';

    % different interest rates
    for ii = [0, 5]
        name = sprintf('r = %g%% p.a.', ii);
        params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
        params(end).r = ii/100;
        params(end).group = {'Q6'};
        params(end).descr = name;
        
        % if ii == 0
        %     params(end).betaH0 = -1e-3;
        %     params(end).beta0 = 0.8;
        % end

        params(end).tex_header = 'r';
        params(end).tex_header_values = {struct('r', sprintf('%g', ii/100))};
    end
    
    % interest rate heterogeneity
    name = 'Permanent r het, r in {0,2,4} p.a.';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).r = [0, 2, 4] / 100;
    params(end).betaH0 = -1e-4;
    params(end).beta0 = 0.973149481985717;
    params(end).group = {'Q6'};
    params(end).descr = 'r in {0, 2, 4}';
    params(end).tex_header = 'r';
    params(end).tex_header_values = {struct('r', '0, 0.02, 0.04')};
    
    name = 'Permanent r het, r in {-2,2,6} p.a.';
    params(end+1) = setup.Params(ifreq,name, quarterly_b_params);
    params(end).r = [-2, 2, 6] / 100;
    params(end).betaH0 = 1e-5;
    params(end).beta0 = 0.960885729527277;
    params(end).group = {'Q6'};
    params(end).descr = 'r in {-2,2,6}';
    params(end).tex_header = 'r';
    params(end).tex_header_values = {struct('r', '-0.02, 0.02, 0.06')};


%         % different tax rates
%         for itax = [0.05, 0.1, 0.15, 0.25]
%             name = [lfreq ' LabTax' num2str(itax)];
%             params(end+1) = setup.Params(ifreq,name,'');
%             params(end).labtaxlow = itax;
%         end


%         % bequest curvature
%         for bcurv = [0.1 0.5 1 2 5]
%             name = [lfreq ' BeqWt0.02 BeqLux0.01 BeqCurv' num2str(bcurv)];
%             params(end+1) = setup.Params(ifreq,name,IncomeProcess);
%             params(end).bequest_weight = 0.02;
%             params(end).bequest_luxury = 0.01;
%             params(end).bequest_curv   = bcurv;target_value
%         end
   
    % fixed beta heterogeneity
    for ibw = [0.001, 0.005, 0.01]
        name = sprintf('Beta5, pSwitch0, pSpacing%g', ibw);
        params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
        params(end).nbeta = 5;
        params(end).betawidth = ibw;
        params(end).prob_zswitch = 0;
        params(end).beta0 = 0.956194383870642;
        params(end).group = {'Q2'};
        params(end).descr = sprintf('p = %g, spacing = %g', 0, ibw);
        params(end).tex_header = '$\beta$ het.';
        params(end).tex_header_values = {struct('pswitch', 0, 'width', ibw)};
        
        if ibw == 0.005
            params(end).betaH0 = -1e-4;
        elseif ibw == 0.01
            params(end).betaH0 = -1e-3;
        end
    end

    % random beta heterogeneity
    for ibw = [0.01]
        for bs = [1/50, 1/10]
            name = sprintf('Beta5, pSwitch%g, pSpacing%g', bs, ibw);
            params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
            params(end).nbeta = 5;
            params(end).betawidth = ibw;
            params(end).prob_zswitch = bs;
            params(end).group = {'Q2'};
            params(end).descr = sprintf('p = %g, spacing = %g', bs, ibw);
            params(end).tex_header = '$\beta$ het.';
            params(end).tex_header_values = {struct('pswitch', bs, 'width', ibw)};

            if bs == 1 / 50
                params(end).betaH0 = -1.2e-2;
            else
                params(end).betaH0 = -1e-3;
            end
        end
    end

    % Different risk aversion coeffs
    for ira = [0.5, 2, 6]
        name = sprintf('CRRA = %g', ira);
        params(end+1) = setup.Params(4, name, quarterly_b_params);
        params(end).risk_aver = ira;
        if (ifreq==4 && ira==4) || ira==6
            params(end).betaL = 0.5;
        end

        if ira == 6
            params(end).betaH0 = -1e-3;
        end
        params(end).group = {'Q3'};
        params(end).descr = 'CRRA';
        params(end).tex_header = 'CRRA';
        params(end).tex_header_values = {struct('riskaver', ira, 'ies', 1 / ira)};
    end
    
    % CRRA with IES heterogeneity
    name = 'CRRA w/IES betw exp(-1), exp(1)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).risk_aver = 1./ exp([-1 -0.5 0 0.5 1]);
    if params(end).freq == 4
        params(end).betaH0 =  - 2e-3;
    end
    params(end).group = {'Q3'};
    params(end).descr = 'RA = exp(1), ..., exp(-1), IES = exp(-1), ..., exp(1)';
    params(end).tex_header = 'CRRA';
    params(end).tex_header_values = {struct('riskaver', 'exp(1), ..., exp(-1)', 'ies', 'exp(-1), ..., exp(1)')};
    
    name = 'CRRA w/IES betw exp(-2), exp(2)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).risk_aver = 1 ./ exp([-2 -1 0 1 2]);
    if params(end).freq == 4
        params(end).betaH0 = -1e-3;
    end
    params(end).beta0 = 0.911905140057402;
    params(end).group = {'Q3'};
    params(end).descr = 'RA = exp(2), ..., exp(-2), IES = exp(-2), ..., exp(2)';
    params(end).tex_header = 'CRRA';
    params(end).tex_header_values = {struct('riskaver', 'exp(2), ..., exp(-2)', 'ies', 'exp(-21), ..., exp(2)')};

    % epstein-zin, quarterly
    ras = [0.5 8  1    1 8];
    ies = [1   1  0.25 2 2];
    for i = 1:5
        ra_i = ras(i);
        ies_i = ies(i);
        name = sprintf('EZ, ra%g, ies%g', ra_i, ies_i);
        params(end+1) = setup.Params(4, name, quarterly_b_params);
        params(end).risk_aver = ra_i;
        params(end).invies = 1 / ies_i;
        params(end).EpsteinZin = true;
        params(end).group = {'Q4'};
        params(end).descr = sprintf('RA = %g, IES = %g', ra_i, ies_i);
        params(end).tex_header = 'EZ';
        params(end).tex_header_values = {struct('riskaver', ra_i, 'ies', ies_i)};

        if i <= 3
            params(end).betaH0 = - 3e-3;
        else
            params(end).betaH0 = - 1e-3;
        end
    end

    % EZ with IES heterogeneity
    name = 'EZ w/ IES betw exp(-1), exp(1)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).invies = 1 ./ exp([-1 -0.5 0 0.5 1]);
    params(end).EpsteinZin = true;
    if (ifreq == 4)
        params(end).betaH0 = - 3e-3;
    end
    params(end).group = {'Q4'};
    params(end).descr = 'RA = 1, IES = exp(-1), ..., exp(1)';
    params(end).tex_header = 'EZ';
    params(end).tex_header_values = {struct('riskaver', 1, 'ies', 'exp(-1), ..., exp(1)')};
    
    name = 'EZ w/ IES betw exp(-2), exp(2)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).invies = 1 ./ exp([-2 -1 0 1 2]);
    params(end).EpsteinZin = true;
    if (ifreq == 4)
        params(end).betaH0 = - 3e-3;
    end
    params(end).group = {'Q4'};
    params(end).descr = 'RA = 1, IES = exp(-2), ..., exp(2)';
    params(end).tex_header = 'EZ';
    params(end).tex_header_values = {struct('riskaver', 1, 'ies', 'exp(-2), ..., exp(2)')};

    % EZ with risk aversion heterogeneity
    name = 'EZ w/ RA betw exp(-2), exp(2)';
    params(end+1) = setup.Params(ifreq, name, quarterly_b_params);
    params(end).invies = 1;
    params(end).risk_aver = exp([-2 -1 0 1 2]);
    params(end).EpsteinZin = 1;
    params(end).betaH0 = - 3e-3;
    params(end).betaL = 0.96;
    params(end).beta0 = 0.99^4;
    params(end).group = {'Q4'};
    params(end).descr = 'RA = exp(-2), ..., exp(2), IES = 1';
    params(end).tex_header = 'EZ';
    params(end).tex_header_values = {struct('riskaver', 'exp(-2), ..., exp(2)', 'ies', 1)};

    % temptation
    for tempt = [0.01 0.05 0.07]
        name = sprintf('Temptation = %g', tempt);
        params(end+1) = setup.Params(4, name, quarterly_b_params);
        params(end).temptation = tempt;
        if tempt == 0.07
            params(end).beta0 = 1;
            params(end).betaH0 = 4e-4;
        elseif  tempt == 0.05
            params(end).betaH0 = -2e-5;
        end
        params(end).group = {'Q5'};
        params(end).tex_header = 'Temptation';
        params(end).tex_header_values = {struct('tempt', sprintf('%g', tempt))};
    end

    name = 'Temptation, uniform in {0, 0.05, 0.1}';
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).temptation = [0 0.05 0.1];
    params(end).beta0 = 0.999083 ^ 4;
    params(end).betaH0 = -4e-5;
    params(end).group = {'Q5'};
    params(end).descr = 'Temptation in {0, 0.05, 0.1}';
    params(end).tex_header = 'Temptation';
    params(end).tex_header_values = {struct('tempt', '0, 0.05, 0.1')};

    %----------------------------------------------------------------------
    % PART 3a, ANNUAL MODEL
    %----------------------------------------------------------------------
    
    % i
    name = 'Annual, no yT';
    params(end+1) = setup.Params(1, name, annual_params);
    params(end).beta0 = 0.99;
    params(end).nyT = 1;
    params(end).sd_logyT = 0;
    params(end).lambdaT = 0;
    params(end).group = {'Q7'};
    params(end).descr = 'No trans shocks';
    params(end).tex_header = 'Annual (i)';
    params(end).tex_header_values = {struct('description', 'no trans shocks')};

    % iv
    name = 'Annual, Carrol';
    params(end+1) = setup.Params(1, name, annual_params);
    params(end).rho_logyP = 0.999;
    params(end).sd_logyP = sqrt(0.015);
    params(end).sd_logyT = sqrt(0.01);
    params(end).group = {'Q7'};
    params(end).descr = 'Carrol process';
    params(end).tex_header = 'Annual (iv)';
    params(end).tex_header_values = {struct('description', 'Carrol process')};
    
%     % v
%     params(end+1) = setup.Params(1,'A a(v) HighPersNotReEst','');
%     params(end).rho_logyP = 0.99;
%     
%     % vi
%     params(end+1) = setup.Params(1,'A a(vi) LowPersNotReEst','');
%     params(end).rho_logyP = 0.9;
% 
%     % vii
%     params(end+1) = setup.Params(1,'A a(vii) HighPersReEst','');
%     params(end).rho_logyP = 0.99;
%     params(end).sd_logyP = sqrt(0.0088);
%     params(end).sd_logyT = sqrt(0.0667);
    
    % viii
    name = 'Annual, high persistence';
    params(end+1) = setup.Params(1, name, shared_params);
    params(end).rho_logyP = 0.995;
    params(end).sd_logyP = sqrt(0.0043);
    params(end).sd_logyT = sqrt(0.0688);
    params(end).lambdaT = 1;
    params(end).group = {'Q7'};
    params(end).descr = 'High persistence';
    params(end).tex_header = 'Annual (viii)';
    params(end).tex_header_values = {struct('description', 'High persistence')};
%     
%     % ix
%     params(end+1) = Params(1,'A a(ix) HighPersNoTransReEst','');
%     params(end).rho_logyP = 0.99;
%     params(end).sd_logyP = sqrt(0.0088);
%     params(end).nyT = 1;
%     params(end).sd_logyT = sqrt(0);
    
    % x
    name = 'Annual, high nyF = 5';
    params(end+1) = setup.Params(1, name, shared_params);
    params(end).rho_logyP = 0.9158;
    params(end).sd_logyP = sqrt(0.0445);
    params(end).sd_logyT = sqrt(0.0479);
    params(end).sd_logyF = sqrt(0.1801);
    params(end).lambdaT = 1;
    params(end).nyF = 5;
    params(end).group = {'Q7'};
    params(end).descr = 'FE heterogeneity';
    params(end).tex_header = 'Annual (x)';
    params(end).tex_header_values = {struct('description', 'FE heterogeneity')};
    
    %----------------------------------------------------------------------
    % PART 3b, QUARTERLY MODEL
    %----------------------------------------------------------------------
    
    % i quarterly_a
    name = 'quarterly_a';
    params(end+1) = setup.Params(4, name, quarterly_a_params);
    params(end).beta0 = 0.984363510593659;
    params(end).group = {'Q8'};
    params(end).descr = 'quart_a';
    params(end).tex_header = 'Quart (i)';
    params(end).tex_header_values = {struct('description', 'quart$_a$')};
    
    % ii
    name = 'KMP';
    params(end+1) = setup.Params(4, name, shared_params);
    params(end).rho_logyP = 0.9879;
    params(end).sd_logyP = sqrt(0.0109);
    params(end).sd_logyT = sqrt(0.0494);
    params(end).group = {'Q8'};
    params(end).descr = 'KMP';
    params(end).tex_header = 'Quart (ii)';
    params(end).tex_header_values = {struct('description', 'KMP')};

    % iii quarterly_c
    name = 'quarterly_c';
    params(end+1) = setup.Params(4, name, quarterly_c_params);
    params(end).beta0 = 0.984363510593659;
    params(end).group = {'Q8'};
    params(end).descr = 'quart_c';
    params(end).tex_header = 'Quart (iii)';
    params(end).tex_header_values = {struct('description', 'quart$_c$')};
    
    % no transitory shocks
    name = 'no trans shocks';
    params(end+1) = setup.Params(4, name, quarterly_b_params);
    params(end).beta0 = 0.984363510593659;
    params(end).nyT = 1;
    params(end).group = {'Q8'};
    params(end).descr = 'no_trans_shocks';
    params(end).tex_header = 'Quart No Trans';
    params(end).tex_header_values = {struct('description', 'quart no trans')};
        
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
        params = setup.Params.select_by_name(params, runopts.name_to_run);
    end

    params.set_run_parameters(runopts);
    params.make_adjustments();

    %----------------------------------------------------------------------
    % ATTACH CALIBRATOR
    %----------------------------------------------------------------------
    if params.calibrating
        calibration = calibrations(params.index);
        calibrator = DHACalibrator(params, calibration.variables,...
            calibration.target_names, calibration.target_values);
        heterogeneity = setup.Prefs_R_Heterogeneity(params);

        if (params.nbeta > 1) && isequal(heterogeneity.ztrans, eye(params.nbeta))
            % Adjust upper bound on beta if using fixed discount factor heterogeneity
            new_betaH = params.betaH - max(heterogeneity.betagrid0);
            params.set("betaH", new_betaH, true);
        end

        beta_bounds = [params.betaL, params.betaH];
        calibrator.set_param_bounds(beta_bounds);
        calibrator.set_handle(params);
        params.set("calibrator", calibrator, true);
    end
end
