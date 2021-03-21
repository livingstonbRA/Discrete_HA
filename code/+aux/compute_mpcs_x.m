function mpcs_cash = compute_mpcs_x(p, basemodel, grdEGP)
    % Computes MPCs over cash-on-hand grid
    mpcs_cash = cell(numel(p.shocks), 1);
    con_base = basemodel.con;
    for ishock = 1:numel(p.shocks)
        shock_size = p.shocks(ishock);

        con_shock = zeros(p.nx, p.nyP, p.nyF, p.nb);
        for ib = 1:p.nb
        for iyF = 1:p.nyF
        for iyP = 1:p.nyP
            cash_shock = grdEGP.x.matrix(:,iyP,iyF,ib) + shock_size;
            con_shock(:,iyP,iyF,ib) = basemodel.coninterp{iyP,iyF,ib}(cash_shock);
        end
        end
        end

        mpcs = (con_shock - con_base) / shock_size;
        mpcs_cash{ishock} = mpcs;
    end
end