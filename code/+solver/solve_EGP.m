function model = solve_EGP(p, grids, heterogeneity,...
    income, futureshock, periods_until_shock, nextmodel, varargin)
    % This function performs the method of endogenous grid points to find
    % saving and consumption policy functions. It also calls 
    % find_stationary() to compute the stationary distribution over states 
    % via direct methods (rather than simulations) and stores the results 
    % in the 'model' structure.
    %
    % To compute the MPCs out of news, it is necessary for the policy function
    % to reflect the expectation of a future shock. For these cases,
    % the policy functions in 'nextmodel' are used. The variable 'nextmpcshock'
    % is nonzero when a shock is expected next period.
    %
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    parser = inputParser;
    addParameter(parser, 'quiet', false);
    parse(parser, varargin{:});
    quiet = parser.Results.quiet;

    %% ----------------------------------------------------
    % USEFUL OBJECTS/ARRAYS
    % -----------------------------------------------------
    nextmpcshock = (periods_until_shock == 1) * futureshock;

    ss_dims = [p.nx, p.nyP, p.nyF, p.nb];
    R_bc = heterogeneity.R_broadcast;

    % If expected future shock is negative, need to raise today's
    % borrowing limit to satisfy future budget constraints with prob one
    % (in all other cases, it should hold that adj_borr_lims = p.borrow_lim).
    tmp = p.borrow_lim - futureshock - income.minnety;
    adj_borr_lims = max(tmp ./ reshape(p.R, 1, []), p.borrow_lim);

    adj_borr_lims_bc = adj_borr_lims(:);
    if numel(adj_borr_lims_bc) == 1
        adj_borr_lims_bc = repmat(adj_borr_lims_bc, p.nb, 1);
    end
    adj_borr_lims_bc = shiftdim(adj_borr_lims_bc, -3);

    svecs_bc = grids.s.vec + (adj_borr_lims_bc - p.borrow_lim);
    svecs_bc_tax = p.compute_savtax(svecs_bc);

    xmat = grids.x.matrix + R_bc .* (adj_borr_lims_bc - p.borrow_lim);

    tempt_expr = heterogeneity.temptation_broadcast;
    tempt_expr = tempt_expr ./ (1 + tempt_expr);

    % Find xprime as a function of s
    xprime_s = R_bc .* svecs_bc + income.netymatEGP + nextmpcshock;

    % Extend to full state space
    newdims = [ss_dims p.nyT];
    [n1, n2, n3, n4, n5] = size(xprime_s);
    dims0 = [n1, n2, n3, n4, n5];
    newdims(newdims == dims0) = 1;
    xprime_s = repmat(xprime_s, newdims);

    %% ----------------------------------------------------
    % CONSTRUCT EXPECTATIONS MATRIX, ETC...
    % -----------------------------------------------------
    Emat = kron(income.ytrans_live, speye(p.nx));

    % Initial guess for consumption function
    tempt_adj = 0.5 * (max(p.temptation) > 0.05);
    con = (max(heterogeneity.r_broadcast, 0.001) + tempt_adj) .* xmat;
    con = con(:);
    con(con<=0) = min(con(con>0));
    con = reshape(con, ss_dims);

    %% ----------------------------------------------------
    % EGP ITERATION
    % ----------------------------------------------------- 
    iter = 1;
    cdiff = 1;
    conlast = con;
    while iter<p.max_iter && cdiff>p.tol_iter
        if iter > 1
            conlast = conupdate;
        end
        iter = iter + 1;

        % interpolate to get c(x') using c(x)

        % c(x')
        c_xp = get_c_xprime(p, grids, xprime_s, ...
            nextmodel, conlast, nextmpcshock, xmat);

        % MUC in current period, from Euler equation
        muc_s = get_marginal_util_cons(...
            p, income, grids, c_xp, xprime_s, R_bc,...
            Emat, heterogeneity.betagrid_broadcast,...
            heterogeneity.risk_aver_broadcast,...
            tempt_expr, svecs_bc);
     
        % c(s)
        con_s = aux.u1inv(heterogeneity.risk_aver_broadcast, muc_s);
        
        % x(s) = s + stax + c(s)
        x_s = svecs_bc + svecs_bc_tax + con_s;

        % interpolate from x(s) to get s(x)
        sav = get_saving_policy(p, grids, x_s, nextmpcshock,...
            R_bc, svecs_bc, xmat);
        sav_tax = p.compute_savtax(sav);

        % updated consumption function, column vec length of
        conupdate = xmat - sav - sav_tax;

        cdiff = max(abs(conupdate(:)-conlast(:)));
        if (mod(iter, 50) == 0) && ~quiet
            disp(['  EGP Iteration ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end
    end
    
    if cdiff>p.tol_iter
        % EGP did not converge, don't find stationary distribution
        AYdiff = 100000;
        model.EGP_cdiff = cdiff;
        return
    end

    model.sav = sav;
    model.con = conupdate;
    model.EGP_cdiff = cdiff;

    % create interpolants from optimal policy functions
    % and find saving values associated with xvals
    % max_sav = (p.borrow_lim - min(income.netymat(:)) - nextmpcshock) ./ p.R;

    model.savinterp = cell(p.nyP,p.nyF,p.nb);
    model.coninterp = cell(p.nyP,p.nyF,p.nb);
    model.coninterp_ext = cell(p.nyP,p.nyF,p.nb);
    model.coninterp_mpc = cell(p.nyP,p.nyF,p.nb);
    for ib = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        model.savinterp{iyP,iyF,ib} = griddedInterpolant(...
            xmat(:,iyP,iyF,ib), model.sav(:,iyP,iyF,ib), 'linear');

        model.coninterp{iyP,iyF,ib} = griddedInterpolant(...
            xmat(:,iyP,iyF,ib), model.con(:,iyP,iyF,ib), 'linear');

        xmin = xmat(1,iyP,iyF,ib);
        cmin = model.con(1,iyP,iyF,ib);
        lbound = 1e-7;

        % Create adjusted interpolants for negative shocks
        model.coninterp_ext{iyP,iyF,ib} = @(x) extend_interp(...
            model.coninterp{iyP,iyF,ib}, x, xmin, cmin, lbound,...
            adj_borr_lims_bc(1,1,1,ib));

        model.coninterp_mpc{iyP,iyF,ib} = @(x) extend_interp(...
            model.coninterp{iyP,iyF,ib}, x, xmin, cmin, -inf,...
            adj_borr_lims_bc(1,1,1,ib));
    end
    end
    end

    model.adj_borr_lims = adj_borr_lims;
    model.adj_borr_lims_bc = adj_borr_lims_bc;
end

function out = extend_interp(old_interpolant, qvals, gridmin,...
    valmin, lb, blim)
    % Applies a correction to the interpolating function since negative
    % shock may push wealth below borrowing limit (only relevant for
    % computing MPCs out of negative shocks)

    out = zeros(size(qvals));
    adj = qvals < gridmin;
    out(~adj) = old_interpolant(qvals(~adj));
    out(adj) = valmin + qvals(adj) - gridmin;
    out(adj) = min(out(adj), qvals(adj)-blim);
    out = max(out, lb);
end

function c_xprime = get_c_xprime(p, grids, xp_s, nextmodel, conlast,...
    nextmpcshock, xmat)
	% find c as a function of x'
	c_xprime = zeros(size(xp_s));
    dims_nx_nyT = [p.nx, 1, 1, 1, p.nyT];

	for ib  = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
    	xp_s_ib_iyF_iyP = xp_s(:,iyP,iyF,ib,:);

        if isempty(nextmodel)
            % Usual method of EGP
            coninterp = griddedInterpolant(...
                xmat(:,iyP,iyF,ib), conlast(:,iyP,iyF,ib), 'linear');
            c_xprime(:,iyP,iyF,ib,:) = reshape(...
                coninterp(xp_s_ib_iyF_iyP(:)), dims_nx_nyT);
        else
            % Solve conditional on an expected future shock using
            % next period's policy functions
            c_xprime(:,iyP,iyF,ib,:) = reshape(...
                nextmodel.coninterp_ext{iyP,iyF,ib}(xp_s_ib_iyF_iyP(:)), dims_nx_nyT);
        end
    end
    end
    end
end

function muc_s = get_marginal_util_cons(...
	p, income, grids, c_xp, xp_s, R_bc,...
    Emat, betagrid_bc, risk_aver_bc, tempt_expr,...
    svecs_bc)

    savtaxrate = (1 + p.savtax .* (svecs_bc >= p.savtaxthresh));

	% First get marginal utility of consumption next period
	muc_c = aux.utility1(risk_aver_bc, c_xp);
    muc_tempt = -tempt_expr .* aux.utility1(risk_aver_bc, xp_s+1e-7);
    mucnext = reshape(muc_c(:) + muc_tempt(:), [], p.nyT);

    % Integrate
    expectation = Emat * mucnext * income.yTdist;
    expectation = reshape(expectation, [p.nx, p.nyP, p.nyF, p.nb]);

    muc_c_today = R_bc .* betagrid_bc .* expectation;
    muc_beq = aux.utility_bequests1(p.bequest_curv, p.bequest_weight,...
        p.bequest_luxury, svecs_bc);

    muc_s = (1 - p.dieprob) * muc_c_today ./ savtaxrate ...
        + p.dieprob * muc_beq;
end

function sav = get_saving_policy(p, grids, x_s, nextmpcshock, R_bc,...
    svecs_bc, xmat)
	% Finds s(x), the saving policy function on the
	% cash-on-hand grid

    sav = zeros(size(x_s));
    xstar = zeros(p.nyP,p.nyF,p.nb);
    for ib  = 1:p.nb
    for iyF = 1:p.nyF
    for iyP = 1:p.nyP
        adj = xmat(:,iyP,iyF,ib) < x_s(1,iyP,iyF,ib);
        sav(adj,iyP,iyF,ib) = svecs_bc(1,1,1,ib);

        savinterp = griddedInterpolant(x_s(:,iyP,iyF,ib),...
            svecs_bc(:,1,1,ib), 'linear');

        sav(~adj,iyP,iyF,ib) = savinterp(...
            xmat(~adj,iyP,iyF,ib));
    end
    end
    end
end