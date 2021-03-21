classdef Prefs_R_Heterogeneity < handle
	% This class stores values, distributions, and
	% transition matrices for heterogeneity in
	% beta, returns, or another variable (z).
	%
	% The z-dimension is set up to accomodate only one source
	% of additional heterogeneity.
	%
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties (SetAccess = private)
		% Grid of betas, centered around zero
		betagrid0;

		% Actual grid of betas = betagrid0 + beta, unless beta_grid_forced
		% is specified in parameters
		betagrid;
		betagrid_broadcast;

		% Length of z-dimension of heterogeneity
		nz;

		% Returns matrices
		R_broadcast;
		r_broadcast;

		% Matrices for risk aversion, inverse IES, and temptation parameter
		risk_aver_broadcast;
		invies_broadcast;
		temptation_broadcast;

		% Distribution and transition matrices for z-hetergeneity
		zdist;
		ztrans;
		zcumdist;
		zcumtrans;
	end

	methods
		function obj = Prefs_R_Heterogeneity(params)
			nEZ = max(numel(params.risk_aver), numel(params.invies));
			dims = [params.nbeta, nEZ, numel(params.r),...
				numel(params.temptation)];

			if sum(dims>1) > 1
				error("Model only allows one source of pref het")
			end
			obj.nz = max(dims);

			obj.initialize_discount_factor(params);
			obj.initialize_IES_heterogeneity(params);
			obj.initialize_temptation_heterogeneity(params);
            obj.initialize_returns_heterogeneity(params);

            obj.create_dists(params);
            obj.zcumdist = cumsum(obj.zdist);
		    obj.zcumtrans = cumsum(obj.ztrans, 2);
		end

		%% -------------------------------------------------------
	    % Discount Factor Heterogeneity
	    % --------------------------------------------------------
		function initialize_discount_factor(obj, params)
		    if isempty(params.beta_grid_forced)
			    bw = params.betawidth;
			    switch params.nbeta
			        case 1
			            obj.betagrid0 = 0;
			        case 2
			            obj.betagrid0 = [-bw/2 bw/2]';
			        case 3
			            obj.betagrid0 = [-bw 0 bw]';
			        case 4
			            obj.betagrid0 = [-3*bw/2 -bw/2 bw/2 3*bw/2]';
			        case 5
			            obj.betagrid0 = [-2*bw -bw 0 bw 2*bw]';
			    end
			    obj.betagrid = params.beta0 + obj.betagrid0;
			else
				% Impose a specific beta grid
				obj.betagrid = params.beta_grid_forced;
			end
			obj.betagrid_broadcast = reshape(...
				obj.betagrid, [1, 1, 1, params.nbeta]);
		end

		%% -------------------------------------------------------
	    % IES Heterogeneity (Epstein-Zin only)
	    % --------------------------------------------------------
	    function initialize_IES_heterogeneity(obj, params)
	    	n_risk_aver = numel(params.risk_aver);
	    	obj.risk_aver_broadcast = reshape(params.risk_aver,...
	    		[1, 1, 1, n_risk_aver]);

	    	n_IES = numel(params.invies);
	    	obj.invies_broadcast = reshape(params.invies,...
	    		[1, 1, 1, n_IES]);
		end

		%% -------------------------------------------------------
	    % Temptation heterogeneity
	    % --------------------------------------------------------
	    function initialize_temptation_heterogeneity(obj, params)
	    	ntempt = numel(params.temptation);
	    	obj.temptation_broadcast = reshape(...
	    		params.temptation, [1, 1, 1, ntempt]);
		end

		%% -------------------------------------------------------
	    % Returns Heterogeneity
	    % --------------------------------------------------------
		function initialize_returns_heterogeneity(obj, params)
            nr = numel(params.r);
			obj.r_broadcast = reshape(params.r, [1 1 1 nr]);
		    obj.R_broadcast = 1 + obj.r_broadcast;
		end

		%% -------------------------------------------------------
	    % Distributions and Transition Matrix
	    % --------------------------------------------------------
	    function create_dists(obj, params)
	    	if ~isempty(params.zdist_forced)
	    		obj.zdist = params.zdist_forced;
	    		switch_prob = 0;
	    	else
	    		switch_prob = params.prob_zswitch / (obj.nz-1);
	    	end

	    	if obj.nz == 1
	    		obj.zdist = 1;
		        obj.ztrans = 1;
		    else
		    	diagonal = (1-params.prob_zswitch) * ones(obj.nz, 1);
		        off_diag = switch_prob * ones(obj.nz);
		        off_diag = off_diag - diag(diag(off_diag));
		        obj.ztrans = off_diag + diag(diagonal);

		        if isempty(params.zdist_forced) && (switch_prob > 0)
					obj.zdist = aux.ergodicdist(obj.ztrans);
				elseif switch_prob == 0
					obj.zdist = ones(obj.nz, 1) / obj.nz;
				end
		    end
	    end
	end
end