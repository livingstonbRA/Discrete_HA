classdef MPCSimulator < handle

	properties (SetAccess = private)
		Nsim;
		Tmax = 4;

		state_rand;
		yPrand;
		zrand;
		yTrand;
		yPindsim;
		yTindsim;
		yFindsim;
		zindsim;
		diesim;

		zcumtrans;

		agrid_long;

		xsim;
		asim;
		csim;
		ssim;
		csim_noshock;

		ygrosssim;
		ynetsim;
		stdev_loggrossy_A;
		stdev_lognety_A;
		mean_grossy_A;

		partitionsize = 1e5;
    	Npartition;

    	adjust_mpc;

    	mpcs = struct();
    	inc_constrained;

    	a1;

        assets_z_rank_corr;
        wpercentiles;
    	wealth_pctile_interpolant;
	end

	methods
		function obj = MPCSimulator(p, distr, heterogeneity)
			obj.Nsim = p.Nsim;
			obj.Npartition = p.Nsim / obj.partitionsize;

			obj.zcumtrans = heterogeneity.zcumtrans;

			obj.state_rand = rand(obj.Nsim,1,'single');
		    obj.yPrand = rand(obj.Nsim,obj.Tmax,'single');
		    obj.zrand = rand(obj.Nsim,obj.Tmax,'single');
		    obj.yTrand = rand(obj.Nsim,obj.Tmax,'single');
		    obj.yPindsim = ones(obj.Nsim,obj.Tmax,'int8');
		    obj.yTindsim = ones(obj.Nsim,obj.Tmax,'int8');
		    obj.yFindsim = ones(obj.Nsim,1,'int8');
		    obj.zindsim = ones(obj.Nsim,1,'int8');

		    dierand = rand(obj.Nsim, obj.Tmax, 'single');
		    obj.diesim = dierand < p.dieprob;

		    wpinterp = griddedInterpolant(...
		        distr.agrid, cumsum(distr.pmf_a), 'pchip', 'nearest');
		    obj.wealth_pctile_interpolant = @(a) 100 * wpinterp(a);
		end

		function simulate(obj, p, income, grids, heterogeneity, basemodel)
			obj.agrid_long = repmat(grids.a.vec, p.nyP*p.nyF*p.nz, 1);
			obj.draw_from_stationary_dist(p, grids, basemodel);
			obj.simulate_exog_transitions(p, income, heterogeneity);
			obj.simulate_decisions(p, grids, basemodel, 0); % baseline
			obj.computeStatistics(p);

			for ishock = 1:6
				obj.simulate_decisions(p, grids, basemodel, ishock);
				obj.computeMPCs(p, ishock);
			end
		end

		function draw_from_stationary_dist(obj, p, grids, basemodel)
			cumdist = cumsum(basemodel.pmf(:));

			% Indexes
		    yPind_trans = repmat(kron((1:p.nyP)', ones(p.nx_DST, 1)), p.nyF*p.nz, 1);
		    yFind_trans = repmat(kron((1:p.nyF)', ones(p.nx_DST*p.nyP, 1)), p.nz, 1);
		    zind_trans = kron((1:p.nz)', ones(p.nx_DST*p.nyP*p.nyF, 1));

			obj.a1 = zeros(obj.Nsim,1);
		    for ip = 1:obj.Npartition
		        partition = obj.partitionsize*(ip-1)+1:obj.partitionsize*ip;

		        % Location of each draw in SSdist
		        [~,ind] = max(obj.state_rand(partition)<=cumdist', [], 2);
		        
		        % (yPgrid,yFgrid,z) indices
		        obj.yPindsim(partition,1) = yPind_trans(ind);
		        obj.yFindsim(partition) = yFind_trans(ind);
		        obj.zindsim(partition,1) = zind_trans(ind);
		        
		        % Initial assets from stationary distribution
		        obj.a1(partition) = obj.agrid_long(ind);
		    end
		end

		function simulate_exog_transitions(obj, p, income, heterogeneity)
			for it = 1:obj.Tmax
		        live = (obj.diesim(:,it) == 0);
		        [~,obj.yTindsim(:,it)] = max(...
		        	obj.yTrand(:,it) <= income.yTcumdist', [], 2);
		        
		        if it > 1
		            if p.ResetIncomeUponDeath == 1
		                [~,obj.yPindsim(live,it)] = max(...
		                	obj.yPrand(live,it) <= income.yPcumtrans(obj.yPindsim(live,it-1),:), [], 2);
		                [~,obj.yPindsim(~live,it)] = max(...
		                	obj.yPrand(~live,it) <= income.yPcumdist', [], 2);
		            else
		                [~,obj.yPindsim(:,it)] = max(...
		                	obj.yPrand(:,it) <= income.yPcumtrans(obj.yPindsim(:,it-1),:), [], 2);
		            end
		            [~,obj.zindsim(:,it)] = max(...
		            	obj.zrand(:,it) <= obj.zcumtrans(obj.zindsim(:,it-1),:), [], 2);
		        end
		    end

		    obj.ygrosssim = income.yPgrid(obj.yPindsim) .*...
            	income.yFgrid(obj.yFindsim) .* income.yTgrid(obj.yTindsim);

		    obj.ynetsim = income.lumptransfer + (1-p.labtaxlow) * obj.ygrosssim...
                        - p.labtaxhigh * max(obj.ygrosssim-income.labtaxthresh, 0);
		end

		function simulate_decisions(obj, p, grids, basemodel, ishock)
			% Location of households that get pushed below bottom of their grid
	        % in first period upon shock
            obj.adjust_mpc = false(obj.Nsim,1);

	        if ishock > 0
	        	shock = p.shocks(ishock);
	        else
	        	shock = 0;
	        end

	        obj.xsim = zeros(obj.Nsim, obj.Tmax);
	        obj.asim = zeros(obj.Nsim, obj.Tmax);
	        obj.csim = zeros(obj.Nsim, obj.Tmax);
	        obj.ssim = zeros(obj.Nsim, obj.Tmax);
			for it = 1:obj.Tmax
	            % Update cash-on-hand          
	            if it == 1
	                obj.xsim(:,1) = obj.a1 + obj.ynetsim(:,1) + shock;
	            else
	                obj.xsim(:,it) = obj.asim(:,it-1) + obj.ynetsim(:,it);
	            end
	            
	            for ib = 1:p.nz
	            for iyF = 1:p.nyF
	            for iyP = 1:p.nyP
	            	idx = obj.yPindsim(:,it)==iyP & obj.yFindsim==iyF & obj.zindsim(:,it)==ib;
	                
	                if (shock < 0) && (it == 1)
	                    below_grid = obj.xsim(:,it) < grids.x.matrix(1,iyP,iyF,ib);
	                    % Bring households pushed below grid back up to grid
	                    idx_below = idx & below_grid;
	                    obj.xsim(idx_below,it) = grids.x.matrix(1,iyP,iyF,ib);
	                    obj.adjust_mpc = obj.adjust_mpc | idx_below;
	                end
	                obj.ssim(idx,it) = basemodel.savinterp{iyP,iyF,ib}(obj.xsim(idx,it));
	            end
	            end
	            end
	            
	            obj.ssim(:,it) = max(obj.ssim(:,it), p.borrow_lim);

	            if it < obj.Tmax
	                if numel(p.r) > 1
	                    obj.asim(:,it) = p.R(obj.zindsim(:,it))' .* obj.ssim(:,it);
	                else
	                    obj.asim(:,it) = p.R * obj.ssim(:,it);
	                end
	                if p.Bequests == 0
	                    % Assets discarded
	                    obj.asim(obj.diesim(:,it)==1,it) = 0;
	                end
	            end
	        end
	        
	        obj.csim = obj.xsim - obj.ssim - p.savtax * max(obj.ssim-p.savtaxthresh, 0);
			
			if ishock == 0
				obj.csim_noshock = obj.csim;
			end
		end

		function computeStatistics(obj, p)
            y_quarter = obj.ynetsim(:,4) * p.freq / 4;
            obj.inc_constrained = struct();
            obj.inc_constrained.a_sixth_Q = mean(obj.asim(:,3) < (y_quarter/6));
            obj.inc_constrained.a_twelfth_Q = mean(obj.asim(:,3) < (y_quarter/12));
            obj.inc_constrained.x_sixth_Q = mean(obj.xsim(:,4) < (y_quarter/6));
            obj.inc_constrained.x_twelfth_Q = mean(obj.xsim(:,4) < (y_quarter/12));

            obj.stdev_loggrossy_A = std(log(sum(obj.ygrosssim(:,1:p.freq),2)));
		    obj.stdev_lognety_A = std(log(sum(obj.ynetsim(:,1:p.freq),2)));
		    obj.mean_grossy_A = mean(sum(obj.ygrosssim(:,1:p.freq),2));
            
            obj.wpercentiles = prctile(obj.asim(:,3), p.percentiles);

		    % Rank-rank correlation between assets and z
		    if p.nz == 1
		    	obj.assets_z_rank_corr = NaN;
		    else
		    	asset_ranks = obj.wealth_pctile_interpolant(obj.asim(:,3));
		    	tmp = corrcoef([double(obj.zindsim(:,3)), asset_ranks]);
                obj.assets_z_rank_corr = tmp(1,2);
		    end
		end

		function computeMPCs(obj, p, ishock)
			shock = p.shocks(ishock);

			% MPC in period 1 out of period 1 shock
			mpcs_1_1 = (obj.csim(:,1) - obj.csim_noshock(:,1)) / shock;
            mpcs_1_1(obj.adjust_mpc,1) = 1;  
            obj.mpcs(ishock).avg_1_1 = mean(mpcs_1_1);

            % MPC in period 2 out of period 1 shock
            mpcs_1_2 = (obj.csim(:,2) - obj.csim_noshock(:,2)) / shock;
            mpcs_1_2(obj.adjust_mpc) = 0;
            obj.mpcs(ishock).avg_1_2 = mean(mpcs_1_2);

            % MPC in period 3 out of period 1 shock
            mpcs_1_3 = (obj.csim(:,3) - obj.csim_noshock(:,3)) / shock;
            mpcs_1_3(obj.adjust_mpc) = 0;
            obj.mpcs(ishock).avg_1_3 = mean(mpcs_1_3);

            % MPC in period 4 out of period 1 shock
            mpcs_1_4 = (obj.csim(:,4) - obj.csim_noshock(:,4)) / shock;
            mpcs_1_4(obj.adjust_mpc) = 0;
            obj.mpcs(ishock).avg_1_4 = mean(mpcs_1_4);

            % Cumulative MPCs over first 4 periods
            mpcs_1_1to4 = mpcs_1_1 + mpcs_1_2 + mpcs_1_3 + mpcs_1_4;
            obj.mpcs(ishock).avg_1_1to4 = mean(mpcs_1_1to4);
		end

		function resultsUpdate = append_results(obj, results)
			resultsUpdate = results;
			resultsUpdate.a_sixth_sim = obj.inc_constrained.a_sixth_Q;
		    resultsUpdate.a_twelfth_sim = obj.inc_constrained.a_twelfth_Q;
		    resultsUpdate.x_sixth_sim = obj.inc_constrained.x_sixth_Q;
		    resultsUpdate.x_twelfth_sim = obj.inc_constrained.x_twelfth_Q;
            resultsUpdate.assets_z_rank_corr = obj.assets_z_rank_corr;

		    resultsUpdate.mpcs_sim = obj.mpcs;
		end
	end
end