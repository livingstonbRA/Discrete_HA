classdef ComparisonDecomp < handle

	properties
		na;
		nthresholds;
		thresholds;

		p0;
		p1;
		stats0;
		stats1;

		pmf0_a;
		pmf1_a;

		RA_mpcs_available;
		mpc0_ra = NaN;
		mpc1_ra = NaN;

		agrid;
		incompatible_params = false;
		same_model_name;

		Empc0_g0;
		Empc0_g1;
		Empc1_g0;
		Empc1_g1;
		Empc1adj_g0 = NaN;

		integral_m0g0_interp;
		integral_m0g1_interp;

		results = struct();
	end

	methods
		function obj = ComparisonDecomp(p0, p1, stats0, stats1)
			obj.p0 = p0;
			obj.p1 = p1;
			obj.stats0 = stats0;
			obj.stats1 = stats1;

			if ~isequal(p0.nx_DST, p1.nx_DST)
				obj.incompatible_params = true;
			elseif ~isequal(p0.abars, p1.abars)
				obj.incompatible_params = true;
			end

			obj.same_model_name = strcmp(p0.name, p1.name);

			obj.RA_mpcs_available = all(...
				~isnan([obj.stats0.mpc_RA.value, obj.stats1.mpc_RA.value]));

			obj.agrid = stats0.agrid;
			obj.na = numel(obj.agrid);

			obj.nthresholds = numel(obj.p0.abars);
			obj.thresholds = obj.p0.abars;

			obj.initialize();
		end

		function initialize(obj)
			import statistics.Statistics.sfill

			nfill = @(varargin) sfill(NaN, varargin{:});
			decomp_name = strcat('Quarterly MPC decomposition of',...
				' ', 'E[MPC] - E[MPC_baseline]');
			obj.results.description = sfill(decomp_name, 'Description');
			obj.results.Em1_less_Em0 = nfill('E[MPC] - E[MPC_baseline]', 3, '$E[m_1]-E[m_0]$');
			obj.results.term1 = nfill('Effect of MPC fcn', 3);
			obj.results.term1a = nfill('Effect of MPC fcn, level', 3);
			obj.results.term1b = nfill('Effect of MPC fcn, shape', 3);
			obj.results.term2 = nfill('Effect of distr', 3);
			obj.results.term3 = nfill('Interaction of MPCs and distr', 3, 'Interaction');

			obj.results.term1_pct = nfill('Effect of MPC fcn (%)', 1, 'Effect of MPC fcn');
			obj.results.term1a_pct = nfill('Effect of MPC fcn (%), level');
			obj.results.term1b_pct = nfill('Effect of MPC fcn (%), shape');
			obj.results.term2_pct = nfill('Effect of distr (%)', 1, 'Effect of distr');
			obj.results.term3_pct = nfill('Interaction of MPCs and distr (%)', 1, 'Interaction');

			obj.results.complete = false;

			for ia = 1:obj.nthresholds
				thresh = obj.thresholds(ia);
				obj.results.term2a(ia) = nfill(...
					sprintf('Effect of distr, HtM (a <= %g)', thresh), 3, '\quad Distr effect, HtM*');
				obj.results.term2b(ia) = nfill(...
					sprintf('Effect of distr, NHtM (a > %g)', thresh), 3, '\quad Distr effect, NHtM*');

				obj.results.term2a_pct(ia) = nfill(...
					sprintf('Effect of distr (%%), HtM (a <= %g)', thresh), 1, '\quad Distr effect, HtM*');
				obj.results.term2b_pct(ia) = nfill(...
					sprintf('Effect of distr (%%), NHtM (a > %g)', thresh), 1, '\quad Distr effect, NHtM*');
			end
		end

		function perform_decompositions(obj, mpcs0, mpcs1)
			if obj.incompatible_params || obj.same_model_name
				return
			end

			if ~all([obj.p0.MPCs, obj.p1.MPCs])
				return
			end

			obj.make_initial_computations(mpcs0, mpcs1);

			obj.results.Em1_less_Em0.value = obj.Empc1_g1 - obj.Empc0_g0;
			benchmark = obj.results.Em1_less_Em0.value / 100;
			if benchmark == 0
				benchmark = NaN;
			end

			% Term 1: Effect of MPC function
			obj.results.term1.value = obj.Empc1_g0 - obj.Empc0_g0;
			obj.results.term1_pct.value = obj.results.term1.value / benchmark;

			if obj.RA_mpcs_available
				% Term 1a: Effect of MPC function, level
				obj.results.term1a.value = obj.Empc1_g0 - obj.Empc1adj_g0;
				obj.results.term1a_pct.value = obj.results.term1a.value / benchmark;

				% Term 1b: Effect of MPC function, shape
				obj.results.term1b.value = obj.Empc1adj_g0 - obj.Empc0_g0;
				obj.results.term1b_pct.value = obj.results.term1b.value / benchmark;
			end

			% Term 2: Effect of distribution
			obj.results.term2.value = obj.Empc0_g1 - obj.Empc0_g0;
			obj.results.term2_pct.value = obj.results.term2.value / benchmark;

			% Term 3: Interaction
			obj.results.term3.value = (obj.Empc1_g1 - obj.Empc0_g1)...
				- (obj.Empc1_g0 - obj.Empc0_g0);
			obj.results.term3_pct.value = obj.results.term3.value / benchmark;

			for ia = 1:obj.nthresholds
				thresh = obj.thresholds(ia);

				% Term 2a: Dist effect for HtM households
				obj.results.term2a(ia).value = obj.integral_m0g1_interp(thresh)...
					- obj.integral_m0g0_interp(thresh);
				obj.results.term2a_pct(ia).value = obj.results.term2a(ia).value / benchmark;

				% Term 2b: Dist effect for NHtM households
				obj.results.term2b(ia).value =...
					(obj.Empc0_g1 - obj.integral_m0g1_interp(thresh))...
					- (obj.Empc0_g0 - obj.integral_m0g0_interp(thresh));
				obj.results.term2b_pct(ia).value = obj.results.term2b(ia).value / benchmark;
			end

			obj.results.complete = true;
		end

		function make_initial_computations(obj, mpcs0, mpcs1)

			pmf0 = obj.stats0.pmf;
			obj.pmf0_a = obj.stats0.pmf_a;
			mpcs0_a = aux.condl_mpcs(mpcs0, pmf0, obj.pmf0_a);

			if ~isequal(obj.stats0.agrid, obj.stats1.agrid)
				% Need to reconstruct mpc's and pmf for model 1 onto
				% baseline grid

				agrid1_orig = obj.stats1.agrid;
				pmf1_orig = obj.stats1.pmf;
				pmf1_a_orig = obj.stats1.pmf_a;
				mpcs1_a_orig = aux.condl_mpcs(...
					mpcs1, pmf1_orig, pmf1_a_orig);

				% First get cdf interpolant
                cdf1_a_interp = griddedInterpolant(...
                	agrid1_orig, cumsum(pmf1_a_orig), 'pchip', 'nearest');

            	% Next get pmf on the baseline grid
            	obj.pmf1_a = zeros(obj.na, 1);
            	obj.pmf1_a(1) = cdf1_a_interp(obj.agrid(1));
            	for ia = 2:obj.na
            		obj.pmf1_a(ia) = cdf1_a_interp(obj.agrid(ia)) - cdf1_a_interp(obj.agrid(ia-1));
            	end
            	obj.pmf1_a = obj.pmf1_a / sum(obj.pmf1_a);

            	% Next get mpc function on baseline grid
            	mpcs1_a_interp = griddedInterpolant(...
            		agrid1_orig, mpcs1_a_orig, 'pchip', 'nearest');
            	mpcs1_a = mpcs1_a_interp(obj.agrid);
            else
            	pmf1 = obj.stats1.pmf;
				obj.pmf1_a = obj.stats1.pmf_a;
				mpcs1_a = aux.condl_mpcs(mpcs1, pmf1, obj.pmf1_a);
			end

			if obj.RA_mpcs_available
				offset = obj.stats1.mpc_RA.value - obj.stats0.mpc_RA.value;
	            mpcs1_adj = mpcs1_a - offset;
	            obj.Empc1adj_g0 = dot(mpcs1_adj, obj.pmf0_a);
	        end

			obj.Empc0_g0 = dot(mpcs0_a, obj.pmf0_a);
			obj.Empc0_g1 = dot(mpcs0_a, obj.pmf1_a);
			obj.Empc1_g0 = dot(mpcs1_a, obj.pmf0_a);
            obj.Empc1_g1 = dot(mpcs1_a, obj.pmf1_a);

            obj.integral_m0g0_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs0_a, obj.pmf0_a, true);
            obj.integral_m0g1_interp = aux.interpolate_integral(...
            	obj.agrid, mpcs0_a, obj.pmf1_a, true);
		end
	end
end