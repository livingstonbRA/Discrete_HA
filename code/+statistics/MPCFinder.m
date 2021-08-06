classdef MPCFinder < handle
	% This class is used for computing MPCs after solving for
	% the policy functions and stationary distribution.

	% Results are storied in 'mpcs', which is filled with NaNs
	% after instantiation. Next, the solve() method is called
	% and the results can be retrieved.

	% 'basemodel' is a structure containing the policy functions
	% and the stationary distribution of the model

	% 'models' is a cell array indexed by (1) the shock index,
	% and (2) the shock period less the current period, plus one. This
	% array contains the policy functions based on expectation
	% of a future shock.

	properties (SetAccess = private)
		Nstates; % total number of states
		r_mat; % matrix of interest rates, if more than 1

		p;
		income;
		grids;

		con_baseline; % baseline consumption
		con_baseline_yT; % baseline consumption before expectation wrt yT
		basemodel; % baseline model (policy fns, etc)
		models; %  models for other shock timings

		xgrid_yT; % cash-on-hand, includes a yT dimension

		fspace; % object used to interpolate back onto agrid

		mpcs = struct(); % results
		loan = struct();
		loss_in_2_years = struct();
	end

	methods
		function obj = MPCFinder(p, income, grids, heterogeneity,...
			basemodel, models)
			import statistics.Statistics.sfill

			obj.Nstates = p.nx_DST * p.nyP * p.nyF * p.nz;
			obj.basemodel = basemodel;
			obj.models = models;
			obj.fspace = fundef({'spli',grids.a.vec,0,1});
			obj.p = p;
			obj.income = income;

			obj.xgrid_yT = repmat(grids.a.vec + income.netymat_broadcast,...
				[1, 1, 1, p.nz, 1]);
			obj.grids = grids;
			obj.r_mat = heterogeneity.r_broadcast;

		    for ishock = 1:6
		    	shock_label = p.shocks_labels{ishock};
		    	shock_tex = split(p.shocks_labels{ishock}, "$");
		    	shock_label_tex = sprintf('%s\\$%s', shock_tex{1}, shock_tex{2});

		    	% mean mpc in period t out of period s shock
		    	obj.mpcs(ishock).quarterly = sfill(NaN,...
		    		sprintf('Quarterly MPC (%%), out of %s', shock_label), 1,...
		    		sprintf('Quarterly MPC (\\%%), out of %s', shock_label_tex));
		    	obj.mpcs(ishock).annual = sfill(NaN,...
		    		sprintf('Annual MPC (%%), out of %s', shock_label), 1,...
		    		sprintf('Annual MPC (\\%%), out of %s', shock_label_tex));

		    	obj.mpcs(ishock).oneperiod = sfill(NaN,...
		    		'MPC, quarterly or annual (\%)', 1, 'MPC, quarterly or annual (\%)');

		    	obj.mpcs(ishock).shock = sfill(shock_label,...
					'Shock description');
				obj.mpcs(ishock).shock_normalized = sfill(p.shocks(ishock),...
					'Shock size as fraction of mean ann inc');
				obj.mpcs(ishock).median_mpc = sfill(NaN,...
					sprintf('Median one-period MPC (%%), out of %s',...
						shock_label));

                obj.mpcs(ishock).quarterly_htm_biweekly = sfill(NaN,...
					sprintf('Quarterly HtM1 MPC (%%), out of %s', shock_label), 1,...
		    		sprintf('Quarterly MPC\textsuperscript{$\\dagger$} (\\%%), out of %s', shock_label_tex));

                obj.mpcs(ishock).quarterly_htm_a_lt_1000 = sfill(NaN,...
					sprintf('Quarterly HtM1 MPC (%%), out of %s', shock_label), 1,...
		    		sprintf('Quarterly MPC\textsuperscript{$\\dagger$} (\\%%), out of %s', shock_label_tex));

                obj.mpcs(ishock).avg_s_t = cell(5, 9);
                for iss = 1:5
                	for itt = 1:9
                		s_t_label = sprintf('Period %d MPC out of %s period %d shock',...
                			itt, shock_label_tex, iss);
                		obj.mpcs(ishock).avg_s_t{iss,itt} = sfill(NaN,...
                			s_t_label, 1, s_t_label);
                	end
                end

		    	% obj.mpcs(ishock).avg_s_t = NaN(5,9);
		    	% mpcs over states for shock in period 1
		    	obj.mpcs(ishock).mpcs_1_t = cell(1,4);
                % fraction of responders
                obj.mpcs(ishock).mpc_pos = NaN(5,5);
                obj.mpcs(ishock).mpc_neg = NaN(5,5);
                obj.mpcs(ishock).mpc0 = NaN(5,5);
                % conditional mean
                obj.mpcs(ishock).avg_s_t_condl = NaN(5,9);
                % median
                obj.mpcs(ishock).median = NaN(5,5);

                obj.mpcs(ishock).avg_1_1to4 = NaN;
                obj.mpcs(ishock).avg_5_1to4 = NaN;

		    	obj.loan.avg = NaN;
		    	obj.loan.mpc_condl = NaN;
		    	obj.loan.mpc_neg = NaN;
		    	obj.loan.mpc0 = NaN;
		    	obj.loan.mpc_pos = NaN;
		    	obj.loan.median = NaN;
		    	obj.loss_in_2_years.avg = NaN;
		    	obj.loss_in_2_years.mpc_condl = NaN;
		    	obj.loss_in_2_years.mpc_neg = NaN;
		    	obj.loss_in_2_years.mpc0 = NaN;
		    	obj.loss_in_2_years.mpc_pos = NaN;
		    	obj.loss_in_2_years.median = NaN;
	    	end
		end

		function solve(obj, p, grids)
			% This function calls the appropriate methods to find
			% the MPCs.

			fprintf('    Computing baseline consumption...\n')
			obj.get_baseline_consumption();

			for ishock = 1:6
				shock_size = p.shocks(ishock);
				if (p.MPCs_news == 0) || (p.EpsteinZin == 1)
					shockperiods = 1;
				else
					shockperiods = [1 2 5];
				end

				for shockperiod = shockperiods
					immediate_shock = 0;
					fprintf('    Computing MPCs out of period %d shock of size %f...\n',...
						shockperiod, shock_size)
					obj.computeMPCs(ishock, shockperiod, immediate_shock);
				end
            end

            if p.MPCs_loan_and_loss == 1
                % $500 loss in 2 years
                immediate_shock = 0;
                fprintf('    Computing MPCs out of anticipated loss in 2 years...\n')
                obj.computeMPCs(1, 9, immediate_shock);

                % $5000 loan for one year
                fprintf('    Computing MPCs out of loan...\n')
                immediate_shock = p.shocks(6);
                obj.computeMPCs(3, 5, immediate_shock);
            end

			obj.compute_cumulative_mpcs();
		end

		function get_baseline_consumption(obj)
			% Each (a,yP,yF) is associated with nyT possible x values, create this
		    % grid here
		    obj.con_baseline_yT = obj.get_policy(obj.xgrid_yT, obj.basemodel);

		    % take expectation wrt yT
		    obj.con_baseline = reshape(obj.con_baseline_yT, [], obj.p.nyT) * obj.income.yTdist;
		end

		function con = get_policy(obj, x_mpc, model)
			% Computes consumption policy function after taking expectation to get
		    % rid of yT dependence
		    con = zeros(size(x_mpc));
		    for ib = 1:obj.p.nz
		    for iyF = 1:obj.p.nyF
		    for iyP = 1:obj.p.nyP
		        x_iyP_iyF_iyT = x_mpc(:,iyP,iyF,ib,:);
		        con_iyP_iyF_iyT = model.coninterp_mpc{iyP,iyF,ib}(x_iyP_iyF_iyT(:));
		        con(:,iyP,iyF,ib,:) = reshape(con_iyP_iyF_iyT,...
		        	[obj.p.nx_DST 1 1 1 obj.p.nyT]);
		    end
		    end
		    end
		end

		function computeMPCs(obj, ishock, shockperiod, immediate_shock)
			shock = obj.p.shocks(ishock);

			for it = 1:shockperiod
				% transition matrix from period 1 to period t
				if it == 1
					trans_1_t = speye(obj.Nstates);
				else
					previous_period = it - 1;
					one_period_transition = obj.transition_matrix_given_t_s(...
						obj.p, shockperiod, previous_period, ishock);
 					trans_1_t = trans_1_t * one_period_transition;
				end

				% cash-on-hand
				transfer_it = (it == shockperiod) * shock ...
					+ (it == 1) * immediate_shock;
				x_mpc = obj.xgrid_yT + transfer_it;

	            % consumption choice given the shock
	            imodel = shockperiod - it + 1;
	            con = obj.get_policy(x_mpc, obj.models{ishock,imodel});

	            % expectation over yT
	            con = reshape(con, [], obj.p.nyT) * obj.income.yTdist;

	            % now compute IMPC(s,t)
	            Econ = trans_1_t * con;
	            Econ_base = obj.basemodel.statetrans^(it-1) * obj.con_baseline;

	            mpcshock = immediate_shock + (immediate_shock == 0) * shock;
	            mpcs = (Econ - Econ_base) / mpcshock;

	            loc_pos = mpcs(:) > 0;
                dist_vec = obj.basemodel.pmf(:);

                mpc_pos = sum(dist_vec(loc_pos));
                mpc_neg = sum(dist_vec(mpcs(:)<0));
                mpc0 = sum(dist_vec(mpcs(:)==0));

                pmf = obj.basemodel.pmf;
                mpc_ss = reshape(mpcs, size(pmf));
                mpcs_sorted = sortrows([mpc_ss(:), pmf(:)]);

                [cdf_m, iu] = unique(cumsum(mpcs_sorted(:,2)), 'last');
                mpc_pct_interp = griddedInterpolant(cdf_m, mpcs_sorted(iu,1),...
                	'pchip', 'nearest');
                mpc_median = mpc_pct_interp(0.5);
        
                if sum(dist_vec(loc_pos)) > 0
                	mpc_condl = dist_vec(loc_pos)' * mpcs(loc_pos) / sum(dist_vec(loc_pos));
                else
                	mpc_condl = NaN;
                end

                mpcsvec = mpcs(:);

           		% Conditional on HtM (own biweekly inc)
       %     		qinc = obj.income.netymat_broadcast  * (obj.p.freq / 4);
       %     		htm_biweekly =  (obj.grids.a.vec ./ qinc) <= (1 / 6);
           		
       %     		[n1, n2, n3, n4] = size(obj.basemodel.pmf);
       %     		dims1 = [n1, n2, n3, n4];
			    % [n1, n2, n3, n4] = size(htm_biweekly);
			    % dims0 = [n1, n2, n3, n4];
			    % dims1(dims1 == dims0) = 1;
			    % htm_biweekly = repmat(htm_biweekly, dims1);

			    % dist_htm_biweekly = kron(obj.income.yTdist, dist_vec);
       %     		dist_htm_biweekly = dist_htm_biweekly(htm_biweekly(:));
       %     		dist_htm_biweekly = dist_htm_biweekly / sum(dist_htm_biweekly);
       %     		mpc_htm_biweekly = dot(dist_htm_biweekly, mpcsvec(htm_biweekly(:)));

           		% Conditional on HtM (a <= $1000)
           		htm_a_lt_1000 =  obj.grids.a.vec <= 1000 / obj.p.numeraire_in_dollars;
           		dims1 = size(obj.basemodel.pmf);
           		dims1(1) = 1;
           		htm_a_lt_1000 = repmat(htm_a_lt_1000, dims1);
           		dist_htm_a_lt_1000 = dist_vec(htm_a_lt_1000(:));
           		dist_htm_a_lt_1000 = dist_htm_a_lt_1000 / sum(dist_htm_a_lt_1000);
           		mpc_htm_a_lt_1000 = dot(dist_htm_a_lt_1000, mpcsvec(htm_a_lt_1000(:)));

           		% Assign values
	            if immediate_shock > 0
	            	obj.loan.avg = obj.basemodel.pmf(:)' * mpcs(:);
	            	obj.loan.mpc_condl = mpc_condl;
	            	obj.loan.mpc_pos = mpc_pos;
	            	obj.loan.mpc0 = mpc0;
	            	obj.loan.mpc_neg = mpc_neg;
	            	obj.loan.median = mpc_median;
                    return;
	            elseif shockperiod <= 5
	            	obj.mpcs(ishock).avg_s_t{shockperiod,it}.value = obj.basemodel.pmf(:)' * mpcs(:);
                    obj.mpcs(ishock).avg_s_t_condl(shockperiod,it) = mpc_condl;
                    obj.mpcs(ishock).mpc_pos(shockperiod,it) = mpc_pos;
                    obj.mpcs(ishock).mpc_neg(shockperiod,it) = mpc_neg;
                    obj.mpcs(ishock).mpc0(shockperiod,it) = mpc0;
                    obj.mpcs(ishock).median(shockperiod,it) = mpc_median;

                    if (it == 1) && (obj.p.freq == 4)
	                    % obj.mpcs(ishock).quarterly_htm_biweekly.value = mpc_htm_biweekly;
	                    obj.mpcs(ishock).quarterly_htm_a_lt_1000.value = 100 * mpc_htm_a_lt_1000;
	                end
	            elseif shockperiod == 9
	            	obj.loss_in_2_years.avg = obj.basemodel.pmf(:)' * mpcs(:);
	            	obj.loss_in_2_years.mpc_condl = mpc_condl;
	            	obj.loss_in_2_years.mpc_pos = mpc_pos;
	            	obj.loss_in_2_years.mpc0 = mpc0;
	            	obj.loss_in_2_years.mpc_neg = mpc_neg;
	            	obj.loss_in_2_years.median = mpc_median;
                    return;
                end
	            
	            if (shockperiod == 1) && (it >= 1 && it <= 4)
	            	% store is = 1 mpcs for decompositions
	                obj.mpcs(ishock).mpcs_1_t{it} = mpcs;
	            end
			end

			if shockperiod == 1
				obj.computeMPCs_periods_after_shock(ishock,...
					shockperiod, trans_1_t);
			end
		end

		function computeMPCs_periods_after_shock(obj,...
			ishock, shockperiod, trans_1_t)

			% transition probabilities from it = is to it = is + 1
			trans_1_t = trans_1_t * obj.transition_matrix_given_t_s(...
				obj.p, shockperiod, shockperiod, ishock);

	        RHScon = obj.basemodel.statetrans^shockperiod * obj.con_baseline(:);
	        LHScon = obj.con_baseline(:);

	        for it = shockperiod+1:9 % it > shockperiod case, policy fcns stay the same in this region
	            mpcs = (trans_1_t * LHScon - RHScon) / obj.p.shocks(ishock);

	            % state transitions are as usual post-shock
	            obj.mpcs(ishock).avg_s_t{shockperiod,it}.value = obj.basemodel.pmf(:)' * mpcs(:);
	            RHScon = obj.basemodel.statetrans * RHScon;
	            LHScon = obj.basemodel.statetrans * LHScon;
	            
	            if (it >= 1) && (it <= 4)
	                obj.mpcs(ishock).mpcs_1_t{it} = mpcs;
	            end
	        end
		end

		function compute_cumulative_mpcs(obj)
			for ishock = 1:6
				vec_1_1to4 = [cell2mat(obj.mpcs(ishock).avg_s_t(1,1:4)).value];
				obj.mpcs(ishock).avg_1_1to4 = sum(vec_1_1to4);
				vec_5_1to4 = [cell2mat(obj.mpcs(ishock).avg_s_t(5,1:4)).value];
				obj.mpcs(ishock).avg_5_1to4 = sum(vec_5_1to4);
				if obj.p.freq == 4
					obj.mpcs(ishock).quarterly.value = 100 * obj.mpcs(ishock).avg_s_t{1,1}.value;
					obj.mpcs(ishock).annual.value = 100 * obj.mpcs(ishock).avg_1_1to4;
					obj.mpcs(ishock).oneperiod.value = obj.mpcs(ishock).quarterly.value;
				else
					obj.mpcs(ishock).quarterly.value = NaN;
					obj.mpcs(ishock).annual.value = 100 * obj.mpcs(ishock).avg_s_t{1,1}.value;
					obj.mpcs(ishock).oneperiod.value = obj.mpcs(ishock).annual.value;
				end

				obj.mpcs(ishock).median_mpc.value = 100 * obj.mpcs(ishock).median(1,1);
			end
		end

		function transition = transition_matrix_given_t_s(...
			obj, p, is, ii, ishock)
			% Computes the transition matrix between t=ii and 
			% t=ii + 1 given shock in period 'is'

			% shocked cash-on-hand
			shock = (is == ii) * p.shocks(ishock);
			x_mpc = obj.xgrid_yT + shock;

			if (ii == is - 1) && (p.shocks(ishock) < 0)
				adj_borr_lim = (obj.p.borrow_lim - p.shocks(ishock) - obj.income.minnety) ./ p.R;
			else
		        adj_borr_lim = p.borrow_lim;
		    end

			% get saving policy function
			imodel = is - ii + 1;
			sav = zeros(p.nx_DST,p.nyP,p.nyF,p.nz,p.nyT);
			reshape_vec = [p.nx_DST 1 1 1 p.nyT];
			for ib = 1:p.nz
			for iyF = 1:p.nyF
			for iyP = 1:p.nyP
				x_iyP_iyF_iyT = reshape(x_mpc(:,iyP,iyF,1,:),[],1);

				if (shock < 0) && (ii == is)
		        	below_xgrid = x_iyP_iyF_iyT < min(obj.xgrid_yT(1,iyP,iyF,1,:));
		            x_iyP_iyF_iyT(below_xgrid) = min(obj.xgrid_yT(1,iyP,iyF,1,:));
                end
		        
                sav_iyP_iyF_iyT = obj.models{ishock,imodel}.savinterp{iyP,iyF,ib}(x_iyP_iyF_iyT);
		        sav(:,iyP,iyF,ib,:) = reshape(sav_iyP_iyF_iyT, reshape_vec);
			end
			end
			end

			% next period's assets conditional on living
			aprime_live = (1+repmat(obj.r_mat, [1 1 1 1 p.nyT])) .* max(sav, adj_borr_lim);

			% interpolate next period's assets back onto asset grid
			asset_interp = funbas(obj.fspace, aprime_live(:));
		    asset_interp = reshape(asset_interp, obj.Nstates, p.nx_DST*p.nyT)...
		    	* kron(speye(p.nx_DST) ,obj.income.yTdist);
		    if p.Bequests == 1
		        interp_death = asset_interp;
		    else
		        interp_death = sparse(obj.Nstates, p.nx_DST);
		        interp_death(:,obj.grids.i0) = 1;
            end
            
            ytrans_live_long = kron(obj.income.ytrans_live, ones(p.nx_DST,1));
            ytrans_death_long = kron(obj.income.ytrans_death, ones(p.nx_DST,1));

		    % now construct transition matrix
		    transition = sparse(obj.Nstates,obj.Nstates);
		    for col = 1:p.nyP*p.nyF*p.nz
		        newblock_live = ytrans_live_long(:,col) .* asset_interp;
		        newblock_death = ytrans_death_long(:,col) .* interp_death;
		        transition(:,p.nx_DST*(col-1)+1:p.nx_DST*col) = ...
		        	(1-p.dieprob)*newblock_live + p.dieprob*newblock_death;
		    end
		end
	end
end