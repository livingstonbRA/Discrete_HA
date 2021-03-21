classdef Statistics < handle

	properties
		pmf;
		pmf_a;
		cdf_a;
		agrid;

		beta_A;
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
		annual_inc_dollars;

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
		end

		function compute_statistics(obj)
			obj.compute_intro_stats();
			obj.add_params();
			obj.construct_distributions();
			obj.compute_percentiles();
			obj.compute_inequality();
			obj.compute_constrained();
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
		function compute_intro_stats(obj)
			obj.beta_A = sfill(obj.p.beta0 ^ obj.freq,...
				'Beta (annualized)', 3);
		    obj.beta_Q = sfill(obj.p.beta0 ^ (obj.freq/4),...
		    	'Beta (quarterly)', 3);
		    
		    tmp = obj.expectation(obj.grdDST.a.matrix);
		    obj.mean_a = sfill(tmp, 'Mean wealth', 3, 'Mean wealth');

		    xdist = obj.model.xdist(:);

		    mean_s = dot(obj.model.sav_x(:), xdist);
		    obj.mean_s = sfill(mean_s, 'Mean s');

		    tmp = dot(obj.model.sav_x(:)==0, xdist);
		    obj.sav0 = sfill(tmp, 's = 0', 3, '$s = 0$');

		    mean_x = dot(obj.model.xvals(:), xdist);
		    obj.mean_x = sfill(mean_x, 'Mean x');

		    % Income
		    mean_y = dot(obj.model.y_x(:) * obj.freq, xdist);
		    obj.mean_gross_y_annual = sfill(mean_y,...
		    	'Mean gross annual income', 3);

		    obj.std_log_gross_y_annual = sfill(NaN,...
			    'Stdev log gross annual income', 3);
			obj.std_log_net_y_annual = sfill(NaN,...
			    'Stdev log net annual income', 3);

			dollars = sprintf('$%g', obj.p.annual_inc_dollars);
			obj.annual_inc_dollars = sfill(dollars,...
			    'Dollar value of mean gross ann inc (numeraire)');

			% RA mpc
			obj.mpc_RA = sfill(NaN, 'MPC (RA)', 3);
		    if (obj.p.nb == 1) && (~obj.p.EpsteinZin) && isequal(obj.p.temptation, 0) && (obj.p.bequest_weight == 0)
		        tmp = (1-obj.p.dieprob) * obj.p.beta0 * obj.p.R;
		        obj.mpc_RA.value = obj.p.R * tmp ^ (-1 / obj.p.risk_aver) - 1;
		    end
		end

		function add_params(obj)
			obj.params = struct();

			if obj.p.freq == 1
				tmp = 'Annual';
			else
				tmp = 'Quarterly';
			end
			obj.params.freq = sfill(tmp,...
				'Frequency');
		end

		function construct_distributions(obj)
			obj.pmf_a = sum(reshape(obj.pmf, obj.nx, []), 2);
		    obj.cdf_a = cumsum(obj.pmf_a);
		end

		function compute_percentiles(obj)
			w_pct = pct_interp(obj.grdDST.a.vec, obj.cdf_a);

			npct = numel(obj.p.percentiles);
		    obj.wpercentiles = cell(1, npct);
			for ip = 1:npct
				pct_at = obj.p.percentiles(ip);

				tmp_b = sprintf('Wealth, %gth pctile', pct_at);
				obj.wpercentiles{ip} = sfill(...
					w_pct(pct_at/100), tmp_b, 3);
			end
			obj.median_a = sfill(w_pct(0.5), 'Median wealth', 3);
		end

		function compute_inequality(obj)
			% Top liquid wealth shares
			cum_share = cumsum(obj.grdDST.a.vec .* obj.pmf_a);
			cum_share = cum_share / obj.mean_a.value;
            
            [cdf_u, iu] = unique(obj.cdf_a, 'last');
			wshare_interp = griddedInterpolant(cdf_u,...
				cum_share(iu), 'pchip', 'nearest');

			tmp = 1 - wshare_interp(0.9);
			obj.w_top10share = sfill(tmp, 'Wealth, top 10% share', 3, 'Wealth, top 10\% share');

			tmp = 1 - wshare_interp(0.99);
			obj.w_top1share = sfill(tmp, 'Wealth, top 1% share', 3, 'Wealth, top 1\% share');
			
			% Gini coefficient
			tmp = aux.direct_gini(obj.grdDST.a.vec, obj.pmf_a);
			obj.wgini = sfill(tmp, 'Gini coefficient, wealth', 3, 'Wealth, Gini coeff');
		end

		function compute_constrained(obj)
			% Constrained by fraction of mean ann inc
			cinterp = constrained_interp(...
	        	obj.grdDST.a.vec, obj.cdf_a);

		    neps = numel(obj.p.epsilon);
		    obj.constrained = cell(1, neps);
		    obj.constrained_pct = cell(1, neps);
		    for ip = 1:neps
				htm = obj.p.epsilon(ip);

				tmp = cinterp(htm);
				obj.constrained{ip} = sfill(tmp,...
					sprintf('a <= %g', htm), 3, sprintf('$a \\leq %g$', htm));

				obj.constrained_pct{ip} = sfill(tmp,...
					sprintf('a <= %g%% mean ann inc', 100 * htm));
			end

			obj.constrained_dollars = {};
			for ip = 1:numel(obj.p.dollar_thresholds)
				htm = obj.p.dollar_thresholds(ip);
				label = obj.p.dollar_threshold_labels{ip};

				tmp = cinterp(htm);
				obj.constrained_dollars{ip} = sfill(tmp,...
					sprintf('a <= %s', label), 3, sprintf('$a \\leq %s$', "\" + label));
			end

			% Wealth / (quarterly earnings) < epsilon
			a_over_inc = obj.grdDST.a.vec ./ ...
				(obj.income.netymat_broadcast * (obj.p.freq / 4));
		    a_over_inc = repmat(a_over_inc, [1, 1, 1, obj.p.nb, 1]);
		    pmf_AY = obj.pmf(:) * shiftdim(obj.income.yTdist, -1);
		    sorted_mat = sortrows([a_over_inc(:), pmf_AY(:)]);

		    cdf_AY = cumsum(sorted_mat(:,2));
		    vals = sorted_mat(:,1);

		    ay_interp = constrained_interp(vals, cdf_AY);

			obj.a_lt_ysixth = sfill(...
				ay_interp(1/6), 'a_i <= y_i / 6', 3, '$a \leq 1 / 6$ own quarterly inc');
			obj.a_lt_ytwelfth = sfill(...
				ay_interp(1/12), 'a_i <= y_i / 12', 3, '$a \leq 1 / 12$ own quarterly inc');
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

% function [values_out, cdf_out, iu] = unique_sort(values_in, pmf_in, iunique)
% 	sorted_mat = sortrows([values_in(:), pmf_in(:)]);
% 	tmp_cdf = cumsum(sorted_mat(:,2));

% 	if iunique == 1
% 		[values_out, iu] = unique(sorted_mat, 'last');
% 		cdf_out = tmp_cdf(iu);
% 	else
% 		[cdf_out, iu] = unique(tmp_cdf, 'first');
% 		values_out = sorted_mat(iu,1);
% 	end
% end