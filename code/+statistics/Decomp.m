classdef Decomp < handle
	properties
		agrid;
		stats;
        p;

		na;
		nthresholds;

		pmf;
		pmf_a;

		integral_interp;
		integral_norisk_interp;
		cdf_interp;

		Empc;
		Empc_norisk;
		mpc_ra;
		mpcs_a;
		mpcs_norisk;
		mean_a;

		results_RA;
		results_norisk;
	end

	methods
		function obj = Decomp(p, adist, stats)
			obj.p = p;
			obj.agrid = adist.agrid;
			obj.na = numel(obj.agrid);
			obj.mean_a = stats.mean_a.value;
			obj.mpc_ra = stats.mpc_RA.value;

			obj.pmf = adist.pmf;
            obj.pmf_a = adist.pmf_a;

			obj.nthresholds = numel(obj.p.abars);

			obj.initialize();
		end

		function perform_decompositions(obj, mpcs, mpcs_norisk)
			obj.mpcs_norisk = mpcs_norisk;

			obj.make_initial_computations(mpcs);
			obj.decomp_RA();
			obj.decomp_norisk();
		end

		function initialize(obj)
			import statistics.Statistics.sfill
			nfill = @(varargin) sfill(NaN, varargin{:});

			obj.results_norisk = struct();
			obj.results_norisk.completed = false;
			obj.results_norisk.term1 = nfill('RA MPC');
			obj.results_norisk.term1_pct = nfill('RA MPC (%)', 1, 'RA MPC (\%)');
			for ia = 1:obj.nthresholds
				thresh = obj.p.abars(ia);
				obj.results_norisk.term2(ia) = nfill(...
					sprintf('HtM effect (a <= %g)', thresh), 3, 'HtM effect');
				obj.results_norisk.term3(ia) = nfill(...
					sprintf('Non-HtM (a <= %g), constraint effect', thresh), 3, 'Non-HtM, constraint effect');
				obj.results_norisk.term4(ia) = nfill(...
					sprintf('Non-HtM (a <= %g), inc risk effect', thresh), 3, 'Non-HtM, inc risk effect');
			end

			obj.results_RA = struct();
			obj.results_RA.completed = false;
			obj.results_RA.Em1_less_mRA = nfill(...
				'E[MPC] - RA MPC');
			obj.results_RA.term1 = nfill(...
				'Effect of MPC function', 3);
			obj.results_RA.term2 = nfill(...
				'Effect of distribution', 3);
			obj.results_RA.term3 = nfill(...
				'Interaction', 3);
		end

		function make_initial_computations(obj, mpcs)
            obj.mpcs_a = aux.condl_mpcs(mpcs, obj.pmf, obj.pmf_a);

            obj.Empc = dot(obj.mpcs_a, obj.pmf_a);
            obj.Empc_norisk = dot(obj.mpcs_norisk, obj.pmf_a);

            obj.integral_interp = aux.interpolate_integral(...
            	obj.agrid, obj.mpcs_a, obj.pmf_a, true);

            obj.integral_norisk_interp = aux.interpolate_integral(...
            	obj.agrid, obj.mpcs_norisk, obj.pmf_a, true);

            integrand = ones(size(obj.pmf_a));
            obj.cdf_interp = aux.interpolate_integral(...
            	obj.agrid, integrand, obj.pmf_a, true);
		end

		function decomp_RA(obj)
			obj.results_RA.Em1_less_mRA.value = obj.Empc - obj.mpc_ra;

			% Interpolant to compute E[MPC|a=E[a]]
			dsupport = obj.pmf_a > 1e-7;
			mpc_a_interp = griddedInterpolant(...
				obj.agrid(dsupport), obj.mpcs_a(dsupport), 'linear');
			mpc_atmean = mpc_a_interp(obj.mean_a);

			% Term 1: Effect of MPC function
			obj.results_RA.term1.value = mpc_atmean - obj.mpc_ra;

			% Term 2: Effect of distribution
			obj.results_RA.term2.value = 0;

			% Term 3: Interaction
			obj.results_RA.term3.value = (obj.Empc - obj.mpc_ra)...
				- (mpc_atmean - obj.mpc_ra);

			obj.results_RA.completed = true;
		end

		function decomp_norisk(obj)
			% Term 1: RA mpc
			obj.results_norisk.term1.value = obj.mpc_ra;
			obj.results_norisk.term1_pct.value = obj.mpc_ra * 100;

			for ia = 1:obj.nthresholds
				thresh = obj.p.abars(ia);
				% Term 2: HtM effect
				obj.results_norisk.term2(ia).value = ...
					obj.integral_interp(thresh) - obj.mpc_ra * obj.cdf_interp(thresh);

				% Term 3: Constraint effect for non-HtM
				obj.results_norisk.term3(ia).value = ...
					(obj.Empc_norisk - obj.integral_norisk_interp(thresh))...
					- obj.mpc_ra * (1 - obj.cdf_interp(thresh));

				% Term 4: Income risk effect for non-HtM
				obj.results_norisk.term4(ia).value = ...
					(obj.Empc - obj.integral_interp(thresh))...
					- (obj.Empc_norisk - obj.integral_norisk_interp(thresh));
			end

			obj.results_norisk.completed = true;
		end
	end
end