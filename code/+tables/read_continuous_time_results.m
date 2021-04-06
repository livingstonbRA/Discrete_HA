function results = read_continuous_time_results(filepath)
    import statistics.Statistics.sfill

    data = load(filepath);
    stats_orig = data.stats;
    
    stats = struct();
    for ishock = 1:6
        quarterly = struct('value', stats_orig.mpcs(ishock).quarterly.value);
        annual = struct('value', stats_orig.mpcs(ishock).annual.value);
        oneperiod = struct('value', stats_orig.mpcs(ishock).quarterly.value);
        stats.mpcs(ishock) = struct('quarterly', quarterly, 'annual', annual,...
            'oneperiod', oneperiod, 'quarterly_htm_a_lt_1000', struct('value', NaN));
    end
    
    stats.beta_A = struct('value', stats_orig.beta_A.value);
    stats.mean_gross_y_annual = struct('value', stats_orig.mean_gross_y_annual.value);
    stats.std_log_gross_y_annual = struct('value', stats_orig.std_log_gross_y_annual.value);
    stats.std_log_net_y_annual = struct('value', stats_orig.std_log_net_y_annual.value);
    stats.mean_a = struct('value', stats_orig.totw.value);
    stats.sav0 = struct('value', stats_orig.sav0.value);
    stats.constrained = {struct('value', stats_orig.constrained{1}.value)};
    
    stats.constrained_dollars = cell(1, 4);
    for ii = 1:4
        stats.constrained_dollars{ii} = struct('value', stats_orig.constrained_dollars{ii}.value);
    end
    stats.a_lt_ysixth = struct('value', stats_orig.w_lt_ysixth.value);
    stats.a_lt_ytwelfth = struct('value', stats_orig.w_lt_ytwelfth.value);
    
    for ii = 1:8
        stats.wpercentiles{ii} = struct('value', stats_orig.wpercentiles{ii}.value);
    end
    
    stats.w_top10share = struct('value', stats_orig.w_top10share.value);
    stats.w_top1share = struct('value', stats_orig.w_top1share.value);
    stats.wgini = struct('value', stats_orig.wgini.value);
    
    decomp_norisk = struct();
    for ii = 1:3
        decomp_norisk.term1_pct(ii) = struct('value', stats_orig.decomp_norisk(ii).term1_pct.value);
        decomp_norisk.term2(ii) = struct('value', stats_orig.decomp_norisk(ii).term2.value);
        decomp_norisk.term3(ii) = struct('value', stats_orig.decomp_norisk(ii).term3.value);
        decomp_norisk.term4(ii) = struct('value', stats_orig.decomp_norisk(ii).term4.value);
    end
    stats.decomp_norisk = decomp_norisk;
    
    stats.decomp_RA = struct();
    stats.decomp_RA.Em1_less_mRA = struct('value', stats_orig.decomp_RA.mean_mpc_diff.value);
    stats.decomp_RA.term1 = struct('value', stats_orig.decomp_RA.term1.value);
    stats.decomp_RA.term2 = struct('value', stats_orig.decomp_RA.term2.value);
    stats.decomp_RA.term3 = struct('value', stats_orig.decomp_RA.term3.value);
    
    results = struct('stats', stats);
end