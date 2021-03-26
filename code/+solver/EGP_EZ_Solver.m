classdef EGP_EZ_Solver < handle
    % This class finds the policy functions using the
    % method of endogenous grid points for the case of
    % Epstein-Zin utility.
    %
    % Brian Livingston, 2020
	% livingstonb@uchicago.edu
    
	properties (SetAccess = private)
        % parameters and grid
        p;
        grids;
        income;
        heterogeneity;
        
		betagrid;
		betastacked;

        % policy and value functions
		con;
		conupdate;
		V;
		Vupdate;
        sav;

		Emat;
        
        EGP_cdiff = 1e5;

        xp_s;
        
        % policy functions of next period cash-on-hand
		c_xp;
		V_xp;

		Vinterp;

		ezval;

		ss_dims;
    	ss_dims_aug;
		repmat_to_state_space;
		reshape_to_state_space;

		sgrid_repeated;
		sgrid_tax;
	end

	methods
		function obj = EGP_EZ_Solver(p, grids, heterogeneity, income)
            obj.p = p;
            obj.grids = grids;
            obj.income = income;
			obj.betagrid = heterogeneity.betagrid;
			obj.heterogeneity = heterogeneity;

		    % initial guess for consumption function, stacked state combinations
		    r_adj = max(obj.p.r, 0.001);
		    obj.con = r_adj * obj.grids.x.matrix;

		    % initial guess for value function
		    obj.V = obj.con;

		    obj.xp_s = p.R * grids.s.matrix + income.netymat_broadcast;

		    obj.ss_dims = [p.nx, p.nyP, p.nyF, p.nb];
		    obj.ss_dims_aug = [obj.ss_dims p.nyT];

		    obj.repmat_to_state_space = ...
		    	@(arr) extend_dims(arr, obj.ss_dims)
		    obj.reshape_to_state_space = ...
		        @(arr) reshape(arr, obj.ss_dims);

		    obj.create_discount_factor_array();
		    obj.create_other_objects(heterogeneity, income);

		    obj.sgrid_repeated = obj.repmat_to_state_space(grids.s.matrix);
		    obj.sgrid_tax = p.compute_savtax(obj.sgrid_repeated);
		end

		function create_discount_factor_array(obj)
			% discount factor matrix, 
		    obj.betastacked = obj.repmat_to_state_space(...
		    	obj.heterogeneity.betagrid_broadcast);
		    obj.betastacked = sparse(diag(obj.betastacked(:)));
		end

		function create_other_objects(obj, heterogeneity, income)
			% construct xpectations operator (conditional on yT)
		    obj.Emat = kron(income.ytrans_live, speye(obj.p.nx));
		end

		function solve(obj, income)
			iter = 1;
			while (iter < obj.p.max_iter) && (obj.EGP_cdiff > obj.p.tol_iter)
				obj.iterate_once(income);

				obj.EGP_cdiff = max(abs(obj.conupdate(:)-obj.con(:)));
		        if mod(iter,50) ==0
		            disp([' EGP Iteration ' int2str(iter)...
		            		' max con fn diff is ' num2str(obj.EGP_cdiff)]);
		        end

		        obj.con = obj.conupdate;
		        obj.V = obj.Vupdate;

				iter = iter + 1;
			end
		end

		function iterate_once(obj,income)
			obj.con = reshape(obj.con, obj.ss_dims);
			obj.V = reshape(obj.V, obj.ss_dims);

			% store c(x') and V(x')
			obj.update_fns_of_xp();

			% matrix of next period muc, muc(x',yP',yF)
			emuc = obj.get_expected_muc(income);

			% current muc(s)
			muc_s = obj.get_current_muc(income, emuc);

			% c(s) by inverting marginal utility
			con_s = muc_s .^ (-1 ./ obj.heterogeneity.invies_broadcast);

			% x(s)
        	x_s = obj.sgrid_repeated +  obj.sgrid_tax + con_s;

	        % s(x) saving policy function
	        obj.sav = obj.get_sav_x_by_interpolating_x_s(x_s);

	        % update consumption
	        savtax = obj.p.compute_savtax(obj.sav);
	        obj.conupdate = obj.grids.x.matrix - obj.sav - savtax;

	        % compute E[V(x)^(1-riskaver)]^(1/(1-riskaver))
	        obj.update_ezval();

	        obj.update_value_fn();
		end

		function update_fns_of_xp(obj)
			obj.c_xp = zeros(obj.ss_dims_aug);
			obj.V_xp = zeros(obj.ss_dims_aug);
			reshape_nx_nyT = @(arr) reshape(arr, [obj.p.nx, 1, 1, 1, obj.p.nyT]);

			for ib  = 1:obj.p.nb
            for iyF = 1:obj.p.nyF
            for iyP = 1:obj.p.nyP
	            % xp_s_ib_iyF_iyP = xp_s(:,iyP,iyF,ib,:);
	            xp_s_ib_iyF_iyP = obj.xp_s(:,iyP,iyF,1,:);
	            coninterp = griddedInterpolant(obj.grids.x.matrix(:,iyP,iyF,ib),...
	            	obj.con(:,iyP,iyF,ib), 'linear');
	            obj.c_xp(:,iyP,iyF,ib,:) = reshape_nx_nyT(coninterp(xp_s_ib_iyF_iyP(:)));
	            obj.Vinterp{iyP,iyF,ib} = griddedInterpolant(...
	            	obj.grids.x.matrix(:,iyP,iyF,ib), obj.V(:,iyP,iyF,ib), 'linear');
	            obj.V_xp(:,iyP,iyF,ib,:) = reshape_nx_nyT(...
	            	obj.Vinterp{iyP,iyF,ib}(xp_s_ib_iyF_iyP(:)));
            end
            end
            end
		end

		function emuc = get_expected_muc(obj, income)
			% nexts period's muc(x',yP',yF)
	        mucnext = obj.c_xp .^ (-obj.heterogeneity.invies_broadcast) ...
	        	.* obj.V_xp .^ (obj.heterogeneity.invies_broadcast ...
	        		- obj.heterogeneity.risk_aver_broadcast);
	       	EyT_mucnext = reshape(mucnext, [], obj.p.nyT) * income.yTdist;
	        
	        % expected muc
	        savtaxrate  = (1+obj.p.savtax.*(obj.sgrid_repeated(:)>=obj.p.savtaxthresh));
	        mu_cons = obj.p.R * obj.betastacked * obj.Emat * EyT_mucnext ./ savtaxrate;
	        mu_cons = obj.reshape_to_state_space(mu_cons);
	        mu_bequest = aux.utility_bequests1(obj.p.bequest_curv, obj.p.bequest_weight,...
	    		obj.p.bequest_luxury, obj.sgrid_repeated);
	        emuc = (1-obj.p.dieprob) * mu_cons + obj.p.dieprob * mu_bequest;
		end

		function muc_s = get_current_muc(obj, income, emuc)
			V_xp_reshaped = reshape(obj.V_xp, [], obj.p.nyT);
			diff_coeffs = obj.heterogeneity.risk_aver_broadcast...
				- obj.heterogeneity.invies_broadcast;
			one_less_riskaver = 1 - obj.heterogeneity.risk_aver_broadcast;

			tmp = exp(obj.Emat * log(V_xp_reshaped) * income.yTdist);
			ezvalnext_ra_equal1 = obj.reshape_to_state_space(tmp) .^ diff_coeffs;

			tmp = obj.V_xp .^ one_less_riskaver;
			tmp = obj.Emat * reshape(tmp, [], obj.p.nyT) * income.yTdist;
			tmp = obj.reshape_to_state_space(tmp);
			ezvalnext_ra_nequal1 = tmp .^ (diff_coeffs ./ one_less_riskaver);

			ezvalnext = zeros(obj.ss_dims);

			equals_one = (obj.heterogeneity.risk_aver_broadcast == 1);
			equals_one = obj.repmat_to_state_space(equals_one);
			ezvalnext(equals_one) = ezvalnext_ra_equal1(equals_one);
			ezvalnext(~equals_one) = ezvalnext_ra_nequal1(~equals_one);

	        muc_s = emuc .* ezvalnext;
	    end

	    function sav = get_sav_x_by_interpolating_x_s(obj, x_s)
	    	% interpolate from x(s) to get s(x)
	        sav = zeros(obj.ss_dims);
	        for ib  = 1:obj.p.nb
	        for iyF = 1:obj.p.nyF
	        for iyP = 1:obj.p.nyP
	            savinterp = griddedInterpolant(x_s(:,iyP,iyF,ib),...
	            	obj.grids.s.matrix(:,iyP,iyF), 'linear');
	            sav(:,iyP,iyF,ib) = savinterp(obj.grids.x.matrix(:,iyP,iyF,1)); 

	            adj = obj.grids.x.matrix(:,iyP,iyF,1) < x_s(1,iyP,iyF,ib);
		        sav(adj,iyP,iyF,ib) = obj.p.borrow_lim;
	        end
	        end
	        end
	    end

	    function update_ezval(obj)
	    	% interpolate adjusted expected value function on x grid
	        ezval_integrand = zeros(obj.ss_dims_aug);
	        for ib = 1:obj.p.nb
	        for iyF = 1:obj.p.nyF
	        for iyP = 1:obj.p.nyP
	            xp_iyP_iyF_ib = obj.xp_s(:,iyP,iyF,1,:);
	            ezval_integrand(:,iyP,iyF,ib,:) = reshape(...
	            	obj.Vinterp{iyP,iyF,ib}(xp_iyP_iyF_ib(:)),...
	            	[obj.p.nx 1 1 1 obj.p.nyT]);
	        end
	        end
	        end

	        ezval_integrand = ezval_integrand...
	        	.^ (1 - obj.heterogeneity.risk_aver_broadcast);

	        % Take expectation over yT
	        ezval_integrand = reshape(ezval_integrand, [], obj.p.nyT)...
	        	* obj.income.yTdist;
	        % Take expectation over (yP,yF,beta)
	        ezval0 = obj.Emat * ezval_integrand;
	        ezval0 = obj.reshape_to_state_space(ezval0);

	        obj.ezval = zeros(obj.ss_dims);

	        equals_one = (obj.heterogeneity.risk_aver_broadcast == 1);
	        equals_one = obj.repmat_to_state_space(equals_one);

	        ezval_ra_equal1 = exp(ezval0);
	        ezval_ra_nequal1 = ezval0 .^...
	        	(1 ./ (1 - obj.heterogeneity.risk_aver_broadcast));
	        obj.ezval(equals_one) = ezval_ra_equal1(equals_one);
	        obj.ezval(~equals_one) = ezval_ra_nequal1(~equals_one);
	    end

	    function update_value_fn(obj)
	    	beta_arr = obj.heterogeneity.betagrid_broadcast;
	    	invies = obj.heterogeneity.invies_broadcast;
	    	ra = obj.heterogeneity.risk_aver_broadcast;

	    	obj.Vupdate = zeros(obj.ss_dims);

	    	equals_one = (invies == 1);
	    	equals_one = obj.repmat_to_state_space(equals_one);

	    	V_ra_1 = obj.conupdate .^ (1 - beta_arr) .* obj.ezval .^ beta_arr;
	    	V_ra_not1 = (1-beta_arr) .* obj.conupdate .^ (1-invies) ...
	    		+ beta_arr .* obj.ezval .^ (1-invies);
	    	V_ra_not1 = V_ra_not1 .^ (1 ./ (1-invies));

	    	obj.Vupdate(equals_one) = V_ra_1(equals_one);
	    	obj.Vupdate(~equals_one) = V_ra_not1(~equals_one);

	        assert(sum(~isfinite(obj.Vupdate(:)))==0)
	        assert(all(obj.Vupdate(:)>=0))
	    end

	    function model = return_model(obj)
	    	model = struct();
	    	model.sav = obj.sav;
		    model.con = reshape(obj.con,[obj.p.nx obj.p.nyP obj.p.nyF obj.p.nb]);
		    model.V = reshape(obj.V,[obj.p.nx obj.p.nyP obj.p.nyF obj.p.nb]);
		    
		    % create interpolants from optimal policy functions
		    % and find saving values associated with xvals
		    model.savinterp = cell(obj.p.nyP,obj.p.nyF,obj.p.nb);
		    model.coninterp = cell(obj.p.nyP,obj.p.nyF,obj.p.nb);
		    model.coninterp_mpc = cell(obj.p.nyP,obj.p.nyF,obj.p.nb);
		    for ib = 1:obj.p.nb
		    for iyF = 1:obj.p.nyF
		    for iyP = 1:obj.p.nyP
		        model.savinterp{iyP,iyF,ib} = ...
		            griddedInterpolant(obj.grids.x.matrix(:,iyP,iyF,ib),...
		            	model.sav(:,iyP,iyF,ib), 'linear');
		        model.coninterp{iyP,iyF,ib} = ...
		            griddedInterpolant(obj.grids.x.matrix(:,iyP,iyF,ib),...
		            	model.con(:,iyP,iyF,ib),'linear'); 

		        xmin = obj.grids.x.matrix(1,iyP,iyF,ib);
		        cmin = model.con(1,iyP,iyF,ib);
		        lbound = 1e-7;

		        model.coninterp_mpc{iyP,iyF,ib} = @(x) extend_interp(...
		            model.coninterp{iyP,iyF,ib}, x, xmin, cmin, -inf,...
		            obj.p.borrow_lim);
		    end
		    end
            end
            
            model.EGP_cdiff = obj.EGP_cdiff;
	    end
	end

end

function out = extend_interp(old_interpolant, qvals, gridmin,...
    valmin, lb, blim)
    out = zeros(size(qvals));
    adj = qvals < gridmin;
    out(~adj) = old_interpolant(qvals(~adj));
    out(adj) = valmin + qvals(adj) - gridmin;

    out(adj) = min(out(adj), qvals(adj)-blim);
    out = max(out, lb);
end

function out = extend_dims(arr, newdims)
    newdims(newdims == size(arr)) = 1;
    out = repmat(arr, newdims);
end