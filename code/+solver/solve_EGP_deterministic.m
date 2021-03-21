function norisk = solve_EGP_deterministic(p, grids,...
    heterogeneity, varargin)
    % This function uses the method of endogenous grid points to find the
    % policy functions of the deterministic model. Output is in the
    % 'norisk' structure.
    % Brian Livingston, 2020
    % livingstonb@uchicago.edu

    parser = inputParser;
    addParameter(parser, 'quiet', false);
    parse(parser, varargin{:});
    quiet = parser.Results.quiet;

    sgrid_bc = repmat(grids.s.vec, 1, p.nb);
    sgrid_tax = p.compute_savtax(sgrid_bc);
    msavtaxrate = (1 + p.savtax .* (sgrid_bc >= p.savtaxthresh));

    Emat = kron(heterogeneity.ztrans, speye(p.nx));
    r_bc = reshape(p.r, 1, []);
    R_bc = reshape(p.R, 1, []);
    risk_aver_bc = reshape(p.risk_aver, 1, []);
    beta_bc = reshape(heterogeneity.betagrid, 1, []);

    tmp = p.temptation ./ (1 + p.temptation);
    tempt_bc = reshape(tmp, 1, []);

    % initial guess for consumption function
    tempt_adj = 0.5 * (max(p.temptation) > 0.05);
    r_bc_adj = max(r_bc, 0.001);
    con = (r_bc_adj + tempt_adj) .* grids.x.matrix_norisk;
    con = con(:);
    con(con<=0) = min(con(con>0));
    con = reshape(con, [p.nx, p.nb]);

    iter = 0;
    cdiff = 1000;
    while (iter <= p.max_iter) && (cdiff > p.tol_iter)
        iter = iter + 1;
        
        conlast = con;

        muc_next = aux.utility1(risk_aver_bc, conlast);
        tempt_next = -tempt_bc .* aux.utility1(...
            risk_aver_bc, grids.x.matrix_norisk);
        beq_next = aux.utility_bequests1(p.bequest_curv, p.bequest_weight,...
            p.bequest_luxury, sgrid_bc);

        expectation = Emat * (muc_next(:) - tempt_next(:));
        emuc_live = R_bc .* beta_bc .* reshape(expectation, [p.nx, p.nb]);

        muc_today = (1 - p.dieprob) * emuc_live ./ msavtaxrate ...
            + p.dieprob .* beq_next;

        con_today = aux.u1inv(risk_aver_bc, muc_today);
        
        cash1 = con_today + sgrid_bc + sgrid_tax;
        
        sav = zeros(p.nx, p.nb);
        for ib = 1:p.nb
            savinterp = griddedInterpolant(cash1(:,ib), grids.s.vec, 'linear');
            sav(:,ib) = savinterp(grids.x.matrix_norisk(:,ib));

            adj = grids.x.matrix_norisk(:,ib) < cash1(1,ib);
            sav(adj,ib) = p.borrow_lim;
        end

        con = grids.x.matrix_norisk - sav - p.compute_savtax(sav);
        
        cdiff = max(abs(con(:)-conlast(:)));
        if ~quiet && ((mod(iter,500) == 0) || (iter == 1))
            fprintf(' EGP for norisk model, iteration %d, norm = %g\n', iter, cdiff)
        end
    end
    norisk.EGP_cdiff = cdiff;
    norisk.con = con;
    norisk.sav = sav;

    if cdiff > p.tol_iter
        fprintf(' No convergence for norisk model.\n')
        norisk.complete = false;
    else
        for ib = 1:p.nb
            norisk.coninterp{ib} = griddedInterpolant(...
                grids.x.matrix_norisk(:,ib), con(:,ib), 'linear');
            norisk.savinterp{ib} = griddedInterpolant(...
                grids.x.matrix_norisk(:,ib), sav(:,ib), 'linear');
        end
        norisk.complete = true;
    end
end