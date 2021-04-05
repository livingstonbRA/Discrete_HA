classdef Params < handle
    % Usage: params = Params(frequency,name,IncomeProcess)
    % After instantiating a Params object, modify any
    % of the parameters using dot notation. The method
    % 'adjust_if_quarterly' must be called at the end since 
    % it is assumed that the discount factor, returns, etc... 
    % are all entered in annual terms.
    %
    % IncomeProcess is the directory of the income grids,
    % or an empty string to generate the income process
    % in the code.
    %
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    properties
        % Identifiers
        name;
        index;
        group;
        descr;
        tex_header;
        tex_header_values = {};
        
        % Data frequency, 1 (annual) or 4 (quarterly)
        freq;

        % Path for income process (or empty string to generate in code)
        IncomeProcess = '';
        
        % Mean annual income, dollars
        numeraire_in_dollars = 72000;

        % Useful functions to convert between numeraire (mean ann inc) and dollars
        convert_to_dollars;
        convert_from_dollars;
        convert_to_dollars_str;

        % Computation
        max_iter = 1e5; % EGP
        tol_iter = 1.0e-7; % EGP
        Nsim = 1e5; % for optional simulation
        Tsim = 400; % for optional simulation
        Nmpcsim = 2e5; % for optional MPC simulation

        % MPC shock sizes
        shocks_dollars = [-1, -500, -5000, 1, 500, 5000];
        shocks;
        shocks_labels;
        
        % Thresholds (in dollars) theta to compute share households with a < theta
        dollar_thresholds = [1000, 2000, 5000, 10000, 25000]; 
        dollar_threshold_labels;

        % Thresholds (in numeraire) theta to compute share households with a < theta
        epsilon = [0, 0.005, 0.01, 0.02, 0.05, 0.1 0.15]; % fraction of mean ann labor income

        % Percentiles to compute
        percentiles = [10, 25, 50, 75, 90, 95, 99, 99.9]; % in percent
        
        % Wealth thresholds (in numeraire) used for decompositions
        abars = [0, 0.01, 0.05];

        % Cash on hand / savings grid parameters
        nx = 250;
        nx_neg = 0;
        nx_DST = 250;
        nx_neg_DST = 0;
        xmax = 500;
        xgrid_par = 0.1; % 1 for linear, 0 for L-shaped
        xgrid_par_neg = 0.4;
        xgrid_term1wt = 0.01;
        xgrid_term1curv = 0.5;
        borrow_lim = 0;
        nbl_adjustment = 0.99;
        
        % Options
        MakePlots = 0;
        Simulate = 0;
        MPCs = 0;
        MPCs_news = 0;
        MPCs_loan_and_loss = 0;
        DeterministicMPCs = 0;
        savematpath = '';
        SaveOutput = false;

        % Annualized return (later adjusted to quarterly if freq = 4)
        r = 0.02;
        R;

        % Turn on annuities
        annuities = false;
        
        % Annualized death probability (later adjusted to quarterly if freq = 4)
        dieprob = 1 / 50;
        
        % Preferences
        EpsteinZin = 0;
        invies = 2.5; % only relevant for Epstein-Zin
        risk_aver = 1;
    	beta0 = 0.98;
    	temptation = 0;

        % Bounds on beta (used when calibrating)
    	betaL = 0.80; % lower bound
        betaH; % theoretical upper bound
        betaH0 = -1e-3; % adjustment to betaH
        
        % Warm glow bequests: bequest weight = 0 is accidental bequests
        bequest_weight = 0;
        bequest_curv = 1;
        bequest_luxury = 0.01;
        Bequests = true; % 1 for wealth left as bequest, 0 for disappears
        
        % Income
    	nyT	= 11;
        yTContinuous = 0; % applies to simulation only and may not work if set to one
        sd_logyT;
        lambdaT = 1; % arrival rate of transitory shock
        nyP = 11;
        sd_logyP;
        rho_logyP;
        nyF = 1;
        sd_logyF = 0;
        ResetIncomeUponDeath = true;
        
        % Taxation
        labtaxlow = 0; % proportional tax
        labtaxhigh = 0; % additional tax on incomes above threshold
        labtaxthreshpc = 0.99; % percentile of earnings distribution where high tax rate kicks in
        savtax = 0; % tax rate on savings
        savtaxthresh = 0; % in terms of numeraire
        compute_savtax; % useful function to compute tax on saving
        lumptransfer = 0; % 

        % Discount factor heterogeneity
        nbeta = 1;
        beta_dist = 1;
        betawidth = 0; % distance between beta's
        beta_grid_forced = []; % overrides all other beta values if used

        % Length of z-dimension heterogeneity
        nz = 1;

        % Probability of switch in z-dimension
        prob_zswitch = 0;
        zdist_forced;

        % Calibration
        calibrating; % boolean
        calibrate_maxiter = 60;
        calibrate_tol = 1e-6;
        calibrator; % DHACalibrator object
    end

    methods
        function obj = Params(frequency, name, addl_params)
        	% Creates params object
            obj.name = name;
            obj.freq = frequency;
            obj.R = 1 + obj.r;
            
            % Default income processes
            if frequency == 1
                obj.sd_logyT = sqrt(0.0494);
                obj.sd_logyP = sqrt(0.0422);
                obj.rho_logyP =0.9525;
            elseif frequency == 4
                % use quarterly_a if quarterly & no IncomeProcess is given
                obj.sd_logyT = sqrt(0.2087);
                obj.sd_logyP = sqrt(0.01080);
                obj.rho_logyP = 0.9881;
            else
                error('Frequency must be 1 or 4')
            end

            % Override default values with addl_params structure
            if nargin >= 3
                pfields = fields(addl_params)';
                for pfield = pfields
                    obj.(pfield{1}) = addl_params.(pfield{1});
                end
            end
        end

        function obj = set_run_parameters(obj, runopts)
        	% Fast option for debugging
            if runopts.fast
                obj.nx_DST = 10;
                obj.nx = 11;
                obj.Nmpcsim = 1e2;
                obj.nyT = 5;
                obj.nyP = 3;
                obj.Tsim = 100;
            end

            obj.MakePlots = runopts.MakePlots;

            obj.savematpath = runopts.savematpath;
            
            % Options
            obj.calibrating = runopts.calibrate;
            obj.Simulate = runopts.Simulate;
            obj.MPCs = runopts.MPCs;
            obj.MPCs_news = runopts.MPCs_news;
            obj.MPCs_loan_and_loss = runopts.MPCs_loan_and_loss;
            obj.DeterministicMPCs = runopts.DeterministicMPCs;
        end
        
        function obj = set_index(obj)
        	% Reset index values to count from 1 to numel(params)
            ind = num2cell(1:numel(obj));
            [obj.index] = deal(ind{:});
        end

        function obj = set(obj, field, new_val, quiet)
            % Sets the value of a parameter.
            %
            % Inputs
            % ------
            %
            % field : A string containing the parameter
            %   name.
            %
            % new_val : The desired value of the parameter.
            %
            % quiet : An optional argument that, when it
            %   evaluates to true, suppresses printing
            %   to the screen.

            field = char(field);
            if ~isprop(obj, field)
                error("Requested field is not an attribute of Params.");
            end

            obj.(field) = new_val;

            if ~exist('quiet', 'var')
                quiet = false;
            end

            % if ~quiet
            %     disp(strcat(field, sprintf(" has been reset to %.9f", new_val)));
            % end
        end

        function make_adjustments(obj)
            % Makes necessary adjustments after Params object has been constructed
            obj.make_frequency_adjustments();
            obj.make_other_adjustments();

            if ~obj.EpsteinZin
                obj.invies = obj.risk_aver;
            end
        end

        function make_frequency_adjustments(obj)
            % Adjusts relevant parameters such as r to a quarterly frequency,
            obj.R = (1 + obj.r) .^ (1 / obj.freq);
            obj.r = obj.R - 1;

            obj.savtax = obj.savtax / obj.freq;
            obj.Tsim = obj.Tsim * obj.freq; % Increase simulation time if quarterly
            obj.beta0 = obj.beta0^(1 / obj.freq);
            obj.dieprob = 1 - (1 - obj.dieprob) ^ (1 / obj.freq);
            obj.prob_zswitch = 1 - (1 - obj.prob_zswitch) ^ (1 / obj.freq);

            obj.betaL = obj.betaL ^ (1 / obj.freq);

            obj.betaH = 1 ./ (max(obj.R) * (1-obj.dieprob));
            obj.betaH = obj.betaH + obj.betaH0;

            obj.lumptransfer = obj.lumptransfer / obj.freq;
        end

        function make_other_adjustments(obj)
            if obj.annuities
                obj.Bequests = false;
                obj.r = obj.r + obj.dieprob;
                obj.R = 1 + obj.r;
            end

            obj.nbeta = max(obj.nbeta, numel(obj.beta_grid_forced));
            obj.compute_savtax =...
                @(sav) obj.savtax * max(sav - obj.savtaxthresh, 0);
            obj.convert_to_dollars = @(num) num * obj.numeraire_in_dollars;
            obj.convert_from_dollars = @(dollars) dollars / obj.numeraire_in_dollars;
            obj.convert_to_dollars_str = @(num) dollar_representation(...
                obj.convert_to_dollars(num));

            obj.shocks = obj.shocks_dollars / obj.numeraire_in_dollars;

            % Create printable labels for shock sizes
            if isempty(obj.shocks_labels)
                obj.shocks_labels = {};
                for ishock = 1:numel(obj.shocks)
                    obj.shocks_labels{ishock} = dollar_representation(obj.shocks_dollars(ishock));
                end
            end

            if obj.EpsteinZin
                obj.DeterministicMPCs = false;
            end

            % Create printable labels for dollar wealth thresholds
            if isempty(obj.dollar_threshold_labels)
                obj.dollar_threshold_labels = {};
                for ii = 1:numel(obj.dollar_thresholds)
                    obj.dollar_threshold_labels{ii} = dollar_representation(obj.dollar_thresholds(ii));
                end
            end

            % Convert dollar thresholds back to units of numeraire
            for ii = 1:numel(obj.dollar_thresholds)
                obj.dollar_thresholds(ii) = obj.dollar_thresholds(ii) / obj.numeraire_in_dollars;
            end

            if isempty(obj.descr)
                obj.descr = obj.name;
            end
        end
    end

    methods (Static)
        function objs = select_by_name(objs, name_to_run)
            % Selects parameterization by name

            name_found = false;
            if ~isempty(name_to_run)
                % Indices of selected names within params
                matches = ismember({objs.name}, {name_to_run});
                if sum(matches) == 1
                    objs = objs(matches);
                    name_found = true;
                end
            end

            if ~name_found
                error('Parameterization name not found or not unique')
            end
        end
    end
end

function dollar_str = dollar_representation(quantity)
    if quantity < 0
        pref = '-$';
    else
        pref = '$';
    end

    dollar_str = sprintf('%s%g', pref, abs(quantity));
end
