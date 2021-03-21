function ymoments = simulate_income_moments(params, income)
    Nsim = 2e5;
    Tburn = 150;

    rng(200);

    %% Simulate
    yFrand = rand(Nsim, 1);
    [~,yFindsim] = max(bsxfun(@lt,yFrand,income.yFcumdist'),[],2);
    yFsim = income.yFgrid(yFindsim);

    yPrand = rand(Nsim, Tburn);
    for it = 1:Tburn
        if it == 1
            [~,yPindsim] = max(bsxfun(@lt,yPrand(:,it),income.yPcumdist'),[],2);
        else
            [~,yPindsim] = max(bsxfun(@lt,yPrand(:,it),income.yPcumtrans(yPindsim,:)),[],2);
        end
    end

    yTrand = rand(Nsim, 4);
    yPrand = rand(Nsim, 4);
    for it = 1:params.freq
        [~, yTindsim] = max(bsxfun(@lt, yTrand(:,it), income.yTcumdist'), [], 2);
        yTsim(:,it) = income.yTgrid(yTindsim);

        [~, yPindsim] = max(bsxfun(@lt, yPrand(:,it), income.yPcumtrans(yPindsim,:)), [], 2);
        yPsim(:,it) = income.yPgrid(yPindsim);
    end

    period_y = yFsim .* yTsim .* yPsim;
    period_nety = params.lumptransfer + (1-params.labtaxlow) * period_y ...
        - params.labtaxhigh * max(period_y-income.labtaxthresh, 0);

    % Annual income
    gross_y = sum(period_y, 2);
    net_y = sum(period_nety, 2);
    
    %% Moments
    ymoments = struct();
    ymoments.std_log_y = std(log(gross_y));
    ymoments.std_log_nety = std(log(net_y));
end