classdef Statistics < handle

	properties
		pmf;
		pmf_a;
		cdf_a;
		agrid;

		xgrid_variables;

		beta_A;
		beta_A_effective;
		beta_Q;

		mean_a;
		mean_s;
		mean_x;
		median_a;
		sav0;

		params;

		mean_gross_y_annual;
		std_log_gross_y_annual;
		std_log_net_y_annual;
		numeraire_in_dollars;

		wpercentiles;

		w_top10share;
		w_top1share;
		wgini;

		mpcs;
		mpc_RA;
		decomp_RA;
		decomp_norisk;

		constrained;
		constrained_pct;
		constrained_dollars;

		a_lt_ysixth;
		a_lt_ytwelfth;
	end

	properties (Access=protected)
		p;
		income;
		grdDST;
		model;

		nx;
		freq;
	end

	methods
		function obj = Statistics(p, income, grdDST, model)
			obj.p = p;
			obj.income = income;
			obj.grdDST = grdDST;
			obj.model = model;

			obj.nx = p.nx_DST;
			obj.freq = p.freq;

			obj.pmf = model.pmf;
			obj.agrid = grdDST.a.vec;

			obj.set_labels();
		end

		function compute_statistics(obj)
			obj.compute_intro_stats();
			obj.add_params();
			obj.construct_distributions();
			obj.compute_percentiles();
			obj.compute_inequality();
			obj.compute_constrained();
			obj.construct_xgrid_variables();
		end

		function add_mpcs(obj, mpcs_obj)
			obj.mpcs = mpcs_obj.mpcs;
		end

		function add_decomps(obj, decomps)
			obj.decomp_RA = decomps.results_RA;
			obj.decomp_norisk = decomps.results_norisk;
		end
	end

	methods (Access=protected)
		function set_labels(obj)
			% Intro stats
			obj.beta_A = empty_stat('Beta (annualized)', 3);
			obj.beta_Q = empty_stat('Beta (quarterly)', 3);
			obj.beta_A_effective = empty_stat('Effective discount rate', 3);
			obj.mean_a = empty_stat('Mean wealth', 3, 'Mean wealth');
			obj.mean_s = empty_stat('Mean s');
			obj.sav0 = empty_stat('s = 0', 3, '$s = 0$');
			obj.mean_x = empty_stat('Mean x');
			obj.mean_gross_y_annual = empty_stat('Mean gross annual income', 3);
			obj.std_log_gross_y_annual = empty_stat(...
				'Stdev log gross annual income', 3);
			obj.std_log_net_y_annual = empty_stat(...
			    'Stdev log net annual income', 3);
			obj.numeraire_in_dollars = empty_stat(...
				'Dollar value of mean gross ann inc (numeraire)');
			obj.mpc_RA = empty_stat('MPC (RA)', 3);

			% Params
			obj.params.freq = empty_stat('Frequency');

			% Percentiles and shares
			npct = numel(obj.p.percentiles);
		    obj.wpercentiles = cell(1, npct);
			for ip = 1:npct
				wpct_lab = sprintf('Wealth, %gth pctile', obj.p.percentiles(ip));
				obj.wpercentiles{ip} = empty_stat(wpct_lab, 3);
			end
			obj.median_a = empty_stat('Median wealth', 3);

			obj.w_top10share = empty_stat(...
				'Wealth, top 10% share', 3, 'Wealth, top 10\% share');
			obj.w_top1share = empty_stat(...
				'Wealth, top 1% share', 3, 'Wealth, top 1\% share');
			obj.wgini = empty_stat('Gini coefficient, wealth', 3, 'Wealth, Gini coeff');

			% Wealth constrained
			neps = numel(obj.p.epsilon);
		    obj.constrained = cell(1, neps);
		    obj.constrained_pct = cell(1, neps);
		    for ip = 1:neps
				htm = obj.p.epsilon(ip);
				obj.constrained{ip} = empty_stat(...
					sprintf('a <= %g', htm), 3, sprintf('$a \\leq %g$', htm));

				obj.constrained_pct{ip} = empty_stat(...
					sprintf('a <= %g%% mean ann inc', 100 * htm));
			end

			neps = numel(obj.p.dollar_thresholds);
			obj.constrained_dollars = cell(1, neps);
			for ip = 1:neps
				label = obj.p.dollar_threshold_labels{ip};
				obj.constrained_dollars{ip} = empty_stat(...
					sprintf('a <= %s', label), 3, sprintf('$a \\leq %s$', "\" + label));
			end

			% HtM stats
			obj.a_lt_ysixth = empty_stat(...
				'a_i <= y_i / 6', 3, '$a \leq 1 / 6$ own quarterly inc');
			obj.a_lt_ytwelfth = empty_stat(...
				'a_i <= y_i / 12', 3, '$a \leq 1 / 12$ own quarterly inc');
		end

		function compute_intro_stats(obj)
			obj.beta_A.value = obj.p.beta0 ^ obj.freq;
		    obj.beta_Q.value = obj.p.beta0 ^ (obj.freq / 4);
		    obj.beta_A_effective.value = obj.beta_A.value ...
		    	* (1 - obj.p.dieprob) ^ obj.freq;

		    obj.mean_a.value = obj.expectation(obj.grdDST.a.matrix);

		    xdist = obj.model.xdist(:);
		    obj.mean_s.value = dot(obj.model.sav_x(:), xdist);
		    obj.sav0.value = dot(obj.model.sav_x(:)==0, xdist);
		    obj.mean_x.value = dot(obj.model.xvals(:), xdist);

		    obj.mean_gross_y_annual.value = dot(obj.model.y_x(:) * obj.freq, xdist);
			obj.numeraire_in_dollars.value = sprintf('$%g', obj.p.numeraire_in_dollars);

		    if (obj.p.nz == 1) && (~obj.p.EpsteinZin) && isequal(obj.p.temptation, 0) ...
		    	&& (obj.p.bequest_weight == 0)
		        tmp = (1-obj.p.dieprob) * obj.p.beta0 * obj.p.R;
		        obj.mpc_RA.value = obj.p.R * tmp ^ (-1 / obj.p.risk_aver) - 1;
		    end
		end

		function add_params(obj)
			obj.params = struct();

			if obj.p.freq == 1
				obj.params.freq.value = 'Annual';
			else
				obj.params.freq.value = 'Quarterly';
			end
		end

		function construct_distributions(obj)
			obj.pmf_a = sum(reshape(obj.pmf, obj.nx, []), 2);
		    obj.cdf_a = cumsum(obj.pmf_a);
		end

		function compute_percentiles(obj)
			w_pct = pct_interp(obj.grdDST.a.vec, obj.cdf_a);
			for ip = 1:numel(obj.p.percentiles)
				obj.wpercentiles{ip}.value = w_pct(obj.p.percentiles(ip) / 100);
			end
			obj.median_a.value = w_pct(0.5);
		end

		function compute_inequality(obj)
			% Top liquid wealth shares
			cum_share = cumsum(obj.grdDST.a.vec .* obj.pmf_a);
			cum_share = cum_share / obj.mean_a.value;
            
            [cdf_u, iu] = unique(obj.cdf_a, 'last');
			wshare_interp = griddedInterpolant(cdf_u,...
				cum_share(iu), 'pchip', 'nearest');

			obj.w_top10share.value = 1 - wshare_interp(0.9);
			obj.w_top1share.value = 1 - wshare_interp(0.99);
			
			% Gini coefficient
			obj.wgini.value = aux.direct_gini(obj.grdDST.a.vec, obj.pmf_a);
		end

		function compute_constrained(obj)
			% Constrained by fraction of mean ann inc
			cinterp = constrained_interp(...
	        	obj.grdDST.a.vec, obj.cdf_a);

		    for ip = 1:numel(obj.p.epsilon)
				htm = obj.p.epsilon(ip);

				obj.constrained{ip}.value = cinterp(htm);
				obj.constrained_pct{ip}.value = cinterp(htm);
			end

			for ip = 1:numel(obj.p.dollar_thresholds)
				htm = obj.p.dollar_thresholds(ip);
				obj.constrained_dollars{ip}.value = cinterp(htm);
			end

			% Wealth / (quarterly earnings) < epsilon
			a_over_inc = obj.grdDST.a.vec ./ ...
				(obj.income.netymat_broadcast * (obj.p.freq / 4));
		    a_over_inc = repmat(a_over_inc, [1, 1, 1, obj.p.nz, 1]);
		    pmf_AY = obj.pmf(:) * shiftdim(obj.income.yTdist, -1);
		    sorted_mat = sortrows([a_over_inc(:), pmf_AY(:)]);

		    cdf_AY = cumsum(sorted_mat(:,2));
		    vals = sorted_mat(:,1);

		    ay_interp = constrained_interp(vals, cdf_AY);

			obj.a_lt_ysixth.value = ay_interp(1/6);
			obj.a_lt_ytwelfth.value = ay_interp(1/12);
		end

		function construct_xgrid_variables(obj)
			obj.xgrid_variables = struct();
            
            dims = [obj.nx*obj.p.nyT, obj.p.nyP, obj.p.nyF, obj.p.nz];
            vars = {'xvals', 'y_x', 'nety_x', 'sav_x', 'con_x', 'xdist'};
            for iv = 1:numel(vars)
                obj.xgrid_variables.(vars{iv}) = zeros(dims);
            end

			% Sort by x values
            dim1 = obj.nx*obj.p.nyT;
            for iyP = 1:obj.p.nyP
                for iyF = 1:obj.p.nyF
                    for iz = 1:obj.p.nz
                        [~, ixvals] = sortrows(obj.model.xvals(:,iyP,iyF,iz));
                        for iv = 1:numel(vars)
                            obj.xgrid_variables.(vars{iv})(:,iyP,iyF,iz) = obj.model.(vars{iv})(ixvals,iyP,iyF,iz);
                        end
                    end
                end
            end
		end

		function out = expectation(obj, vals)
			out = dot(obj.pmf(:), vals(:));
		end
	end

	methods (Static)
		function out = sfill(value, label, varargin)
			out = sfill(value, label, varargin{:});
		end
	end
end

function out = empty_stat(varargin)
	out = sfill(NaN, varargin{:});
end

function out = sfill(value, label, decimals, tex_label)
	if (nargin < 3)
		decimals = 1;
	end

	if (nargin < 4)
		tex_label = label;
	end

	out = struct(...
		'value', value,...
		'label', label,...
		'tex_label', tex_label,...
		'decimals', decimals...
	);
end

function interp_out = pct_interp(values, cdf_x)
	[cdf_x_u, iu] = unique(cdf_x(:), 'first');
	values_u = values(iu);

	if numel(cdf_x_u) >= 2
		interp_out = griddedInterpolant(...
			cdf_x_u, values_u, 'pchip', 'nearest');
	else
		interp_out = @(x) NaN;
	end
end

function interp_out = constrained_interp(values, cdf_x)
	[values_u, iu] = unique(values, 'last');
	cdf_x_u = cdf_x(iu);

	if numel(values_u) >= 2
		interp_out = griddedInterpolant(...
			values_u, cdf_x_u, 'pchip', 'nearest');
	else
		interp_out = @(x) NaN;
	end
end