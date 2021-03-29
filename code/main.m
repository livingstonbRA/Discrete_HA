function results = main(p, varargin)
    % Endogenous Grid Points with AR1 + IID Income
    % Cash on Hand as State variable
    % Includes NIT and discount factor heterogeneity
    
    % This is the main function file for this code repository. Given a
    % structure of parameters, p, this script calls functions primarily to 
    % compute policy functions via the method of endogenous grip points, 
    % and to find the implied stationary distribution over the state space.

    results = struct('norisk',[],'sim',[]);

    parser = inputParser;
    addOptional(parser, 'iterating', false);
    parse(parser, varargin{:});
    iterating = parser.Results.iterating;

    %% --------------------------------------------------------------------
    % HETEROGENEITY IN PREFERENCES/RETURNS
    % ---------------------------------------------------------------------
    heterogeneity = setup.Prefs_R_Heterogeneity(p);
    p.set("nb", heterogeneity.nz, true);

    %% --------------------------------------------------------------------
    % INCOME
    % ---------------------------------------------------------------------
    income = setup.Income(p, heterogeneity);

    %% --------------------------------------------------------------------
    % ASSET GRIDS
    % ---------------------------------------------------------------------
    NBL = -min(income.netymat(:)) / max(p.r);
    loose_constraint = p.nbl_adjustment * NBL;
    if p.borrow_lim <= -1e10
        p.set("borrow_lim", loose_constraint, false);
    end

    % grids for method of EGP
    grdEGP = setup.Grid(p, income, 'EGP');

    % grids for finding stationary distribution
    grdDST = setup.Grid(p, income, 'DST');

    %% --------------------------------------------------------------------
    % MODEL SOLUTION
    % ---------------------------------------------------------------------
    % Get policy functions and stationary distribution for final beta, in
    % 'basemodel' structure
    if p.EpsteinZin
        egp_ez_solver = solver.EGP_EZ_Solver(p, grdEGP, heterogeneity, income);
        egp_ez_solver.solve(income);
        basemodel = egp_ez_solver.return_model();
    else
        nextmpcshock = 0;
        periods_until_shock = 0;
        basemodel = solver.solve_EGP(...
            p, grdEGP, heterogeneity, income, nextmpcshock,...
            periods_until_shock, [], 'quiet', iterating);
    end
    basemodel = solver.find_stationary_adist(...
        p, basemodel, income, grdDST, heterogeneity, 'quiet', iterating);

    if basemodel.EGP_cdiff > p.tol_iter
        % EGP did not converge for beta, escape this parameterization
        return
    end
    
    %% --------------------------------------------------------------------
    % STATISTICS
    % ---------------------------------------------------------------------
    results.stats = statistics.Statistics(p, income, grdDST, basemodel);
    results.stats.compute_statistics();

    ymoments = statistics.simulate_income_moments(p, income);
    results.stats.std_log_gross_y_annual.value = ymoments.std_log_y;
    results.stats.std_log_net_y_annual.value = ymoments.std_log_nety;
    
    %% --------------------------------------------------------------------
    % MPCs FOR MODEL WITHOUT INCOME RISK
    % ---------------------------------------------------------------------
    % Fill deterministic model MPCs with NaNs
    results.norisk.mpcs1_a_direct = cell(1,6);
    for im = 1:6
        results.norisk.mpcs1_a_direct{im} = NaN;
    end

    if p.DeterministicMPCs
        % Try to solve deterministic model
        norisk = solver.solve_EGP_deterministic(...
            p, grdEGP, heterogeneity);

        if norisk.complete
            % If converged, compute MPCs
            results.norisk.mpcs1_a_direct = ...
                statistics.direct_MPCs_by_computation_norisk(...
                    p, norisk, income, heterogeneity, grdDST);
        else
            p.set('DeterministicMPCs', false, true);
        end
    end

    %% --------------------------------------------------------------------
    % SIMULATIONS
    % ---------------------------------------------------------------------
    if p.Simulate
        results.sim = solver.simulate(...
            p, income, basemodel, grdDST, heterogeneity);
    end

    %% --------------------------------------------------------------------
    % DIRECTLY COMPUTED MPCs, IMPC(s,t) where s = shock period
    % ---------------------------------------------------------------------
    maxT = 1;
    if ((p.MPCs_news == 1) || (p.MPCs_loan_and_loss == 1)) && ~p.EpsteinZin
        disp('Solving for policy functions of anticipated future shocks')
        if p.freq == 4
            maxT = 10;
        else
            maxT = 5;
        end
    end
    
    % mpcmodels{ishock,tlshock} stores the policy functions associated with the case
    % where the shock occurs in (tlshock - 1) periods after the current period
    mpcmodels = cell(6, maxT);
    for ishock = 1:6
        % policy functions are the same as baseline when shock is received in
        % the current period
        mpcmodels{ishock,1} = basemodel;

        % get consumption functions conditional on future shock
        for tlshock = 2:maxT
            nextmpcshock = (tlshock == 2) * p.shocks(ishock);
            mpcmodels{ishock,tlshock} = solver.solve_EGP(...
                p, grdEGP, heterogeneity, income, nextmpcshock,...
                tlshock-1, mpcmodels{ishock,tlshock-1}, 'quiet', iterating);
        end
    end

    mpc_finder = statistics.MPCFinder(p, income, grdDST, heterogeneity,...
        basemodel, mpcmodels);
    if p.MPCs
        disp('Computing MPCs')
        mpc_finder.solve(p, grdDST);
    end

    results.stats.add_mpcs(mpc_finder);
    results.mpcs = mpc_finder.mpcs;
    results.mpcs_loan = mpc_finder.loan;
    results.mpcs_loss_in_2_years = mpc_finder.loss_in_2_years;
    clear mpc_finder
    
    %% --------------------------------------------------------------------
    % MPCs via DRAWING FROM STATIONARY DISTRIBUTION AND SIMULATING
    % ---------------------------------------------------------------------
    mpc_simulator = statistics.MPCSimulator(p, results.stats, heterogeneity);
    mpc_simulator.simulate(p, income, grdDST, heterogeneity, basemodel);
    results.mpcs_sim = mpc_simulator.mpcs;
    clear mpc_simulator

    %% --------------------------------------------------------------------
    % DECOMPOSITION 1 (DECOMP OF E[mpc])
    % ---------------------------------------------------------------------
    decomp = statistics.Decomp(p, results.stats, results.stats);
    doDecomposition = (p.nb==1) && (~p.EpsteinZin) && (p.MPCs)...
        && (p.bequest_weight==0) && isequal(p.temptation,0) && (numel(p.r)==1)...
        && (p.DeterministicMPCs);

    if doDecomposition
        mpcs_baseline = reshape(results.mpcs(5).mpcs_1_t{1}, p.nx_DST, []);
        mpcs_norisk = reshape(results.norisk.mpcs1_a_direct{5}, p.nx_DST, []);
        decomp.perform_decompositions(mpcs_baseline, mpcs_norisk);
    end

    results.stats.add_decomps(decomp);
    clear decomp
    
    %% --------------------------------------------------------------------
    % FIGURES
    % ---------------------------------------------------------------------
    % if p.MakePlots
    %     % plot(grdDST.a.vec, cumsum(results.direct.agrid_dist))
    %     % xlim([0 0.2])
        
    %     % Wealth at low yP
    %     nbins = 100;
    %     amin = grdDST.a.vec(1);
    %     amax = {1};
    %     amax_visible = 0.5;

    %     iyP = 1:11;
    %     nyP = numel(iyP);

    %     pmf_a = results.direct.adist(:,iyP,:,:);
    %     pmf_a = pmf_a(:) / sum(pmf_a(:));
    %     pmf_a = reshape(pmf_a, [], nyP);
    %     pmf_a = sum(pmf_a, 2);

    %     wealth_plotter = statistics.WealthPlotter(p, grdDST.a.vec, pmf_a);
    %     [ax, wealth_hist] = wealth_plotter.create_histogram(nbins, amax{:});
    %     title("Wealth distribution, truncated above")
    %     ax.XLim = [amin, amax_visible];
    %     ax.YLim = [0, max(wealth_hist.Values(1:end-1))];

    %     figpath = fullfile('output', 'wealth_distribution.jpg');
    %     saveas(gcf, figpath)
        
    %     %% MPCs Function
    %     fontsize = 12;
    %     mpcs = results.direct.mpcs(5).mpcs_1_t{1};
    %     mpc_plotter = statistics.MPCPlotter(p, grdDST.a.matrix, mpcs);
    %     mpc_plotter.fontsize = fontsize;
    %     mpc_plotter.show_grid = 'on';

    %     yP_indices = [3, 8];
    %     zoomed_window = true;
    %     shock_size = 0.01;
    %     [ax_main, ax_window] = mpc_plotter.create_mpcs_plot_yPs(...
    %                 yP_indices, zoomed_window, shock_size);
    %     ylim_main = ax_main.YLim;

    %     imedian = find(p.percentiles == 50);
    %     median_wealth = results.direct.wpercentiles(imedian);
    %     ax_main = mpc_plotter.add_median_wealth(ax_main, median_wealth);

    %     ax_main.XLim = [0, 5];
    %     ax_main.YLim = ylim_main;

    %     window_max_x = 0.3;
    %     ax_window.YLim = ax_main.YLim;
    %     ax_window.XLim = [0, window_max_x];
    %     xticks(ax_window, [0:0.1:window_max_x])
    %     yticks(ax_window, [0:0.1:0.3])
    %     set(ax_window, 'FontSize', fontsize-2)
    %     ax_window.YTick = ax_main.YTick(1:2:end);

    %     figpath = fullfile('output', 'mpc_function_yPs.jpg');
    %     saveas(gcf, figpath)

    %     %% MPCs Function For Diff Shock Sizes
    %     fontsize = 12;
    %     mpcs = {    results.direct.mpcs(2).mpcs_1_t{1}
    %                 results.direct.mpcs(3).mpcs_1_t{1}
    %                 results.direct.mpcs(5).mpcs_1_t{1}
    %                 results.direct.mpcs(6).mpcs_1_t{1}
    %            };
           
    %     for ii = 1:numel(mpcs)
    %         mpcs{ii} = reshape(mpcs{ii}, [p.nx_DST p.nyP p.nyF p.nb]);
    %     end

    %     mpc_plotter = statistics.MPCPlotter(p, grdDST.a.matrix, mpcs);
    %     mpc_plotter.fontsize = fontsize;
    %     mpc_plotter.show_grid = 'on';

    %     iyP = median(1:p.nyP);
    %     ishocks = [2 3 5 6];
    %     zoomed_window = true;
    %     shock_size = 0.01;
    %     [ax_main, ax_window] = mpc_plotter.create_mpc_plot_shocks(...
    %                 iyP, zoomed_window, ishocks);
    %     ylim_main = ax_main.YLim;

    %     imedian = find(p.percentiles == 50);
    %     median_wealth = results.direct.wpercentiles(imedian);
    %     ax_main = mpc_plotter.add_median_wealth(ax_main, median_wealth);

    %     ax_main.XLim = [0, 5];
    %     ax_main.YLim = ylim_main;

    %     window_max_x = 0.3;
    %     ax_window.YLim = ax_main.YLim;
    %     ax_window.XLim = [0, window_max_x];
    %     xticks(ax_window, [0:0.1:window_max_x])
    %     yticks(ax_window, [0:0.1:0.3])
    %     set(ax_window, 'FontSize', fontsize-2)
    %     ax_window.YTick = ax_main.YTick(1:2:end);

    %     figpath = fullfile('output', 'mpc_function_shocks.jpg');
    %     saveas(gcf, figpath)
    % end
    
    % convert Params object to structure for saving
    results.stats = aux.to_structure(results.stats);
    p.set('calibrator', [], false);
    Sparams = aux.to_structure(p);
    converged = iterating;

    if ~iterating
        save(p.savematpath, 'Sparams', 'results', 'converged')
    end
end