function decomp = borrlim_decomposition(p_baseline, results_baseline,...
	p_no_bc, results_no_bc, return_nans)
    % Performs decompositions comparing models with and without borrowing constraints

	% Construct agrid from baseline parameters
	agrid = results_baseline.direct.agrid;

	% Initialize
	decomp = initialize_to_nan(p_baseline);
    if return_nans
        return
    end

	% Check if required MPCs are available
    baseline_mpcs_available = p_baseline.MPCs && p_baseline.DeterministicMPCs;
    no_bc_mpcs_available = p_no_bc.MPCs && p_no_bc.DeterministicMPCs;
    if ~(baseline_mpcs_available && no_bc_mpcs_available)
        return
    end

    % Dimensions other than assets
    n_het = p_baseline.nyP * p_baseline.nyF * p_baseline.nz;

    %% --------------------------------------------------------------------
    % RA model WITHOUT borrowing constraint
    % ---------------------------------------------------------------------
    RA = struct();
    tmp = (1-p_baseline.dieprob) * p_baseline.0 * p_baseline.R;
    RA.mpc = p_baseline.R * tmp ^ (-1/p_baseline.risk_aver) - 1;

    %% --------------------------------------------------------------------
    % RA model WITH borrowing constraint
    % ---------------------------------------------------------------------
    RA_with_BC = struct();

	mpcs = results_baseline.norisk.mpcs1_a_direct{5};
    mpcs = reshape(mpcs, [p_baseline.nx_DST, 1, 1, p_baseline.nz]);
    mpcs = repmat(mpcs, [1, p_baseline.nyP, p_baseline.nyF, 1]);
    RA_with_BC.mpcs = mpcs;

    %% --------------------------------------------------------------------
    % HA model WITH borrowing constraint (baseline)
    % ---------------------------------------------------------------------
    HA_with_BC = struct();

    HA_with_BC.Empc = results_baseline.direct.mpcs(5).avg_s_t(1,1);

    mpcs = results_baseline.direct.mpcs(5).mpcs_1_t{1};
    HA_with_BC.mpcs = reshape(mpcs, [], n_het);
    HA_with_BC.pmf = reshape(results_baseline.direct.adist, [], n_het);

    % Interpolant for cdf
    HA_with_BC.cdf_interp = cell(n_het, 1);
    for ii = 1:n_het
        sortedAssetDistribution = sortrows([agrid HA_with_BC.pmf(:,ii)]);
        [assetVals, uniqueInds] = unique(sortedAssetDistribution(:,1), 'last');
        cumg = cumsum(sortedAssetDistribution(:,2));
        HA_with_BC.cdf_interp{ii} = griddedInterpolant(assetVals, cumg(uniqueInds), 'linear');
    end

    %% --------------------------------------------------------------------
    % HA model WITHOUT borrowing constraint (loose borrowing constraint)
    % ---------------------------------------------------------------------
    HA = struct();

    mpcs = results_no_bc.direct.mpcs(5).mpcs_1_t{1};
    mpcs = reshape(mpcs, [], n_het);
    ind0 = p_no_bc.nx_neg_DST + 1;
    HA.mpcs = mpcs(ind0:end,:);

    %% --------------------------------------------------------------------
    % Prepare interpolants for integrals over MPCs from zero to eps
    % ---------------------------------------------------------------------
    % For integral over mpcs for HA model with budget constraint
    HA_with_BC.mpc_integral_interp = cell(n_het, 1);
    for ii = 1:n_het
        HA_with_BC.mpc_integral_interp{ii} = aux.interpolate_integral(...
            agrid, HA_with_BC.mpcs(:,ii), HA_with_BC.pmf(:,ii));
    end

    % For integral over mpcs for HA model without budget constraint
    HA.mpc_integral_interp = cell(n_het, 1);
    for ii = 1:n_het
        HA.mpc_integral_interp{ii} = aux.interpolate_integral(...
            agrid, HA.mpcs(:,ii), HA_with_BC.pmf(:,ii));
    end

    % For integral over mpcs for RA model with budget constraint
    RA_with_BC.mpc_integral_interp = cell(n_het, 1);
    for ii = 1:n_het
        RA_with_BC.mpc_integral_interp{ii} = aux.interpolate_integral(...
            agrid, RA_with_BC.mpcs(:,ii), HA_with_BC.pmf(:,ii));
    end

    %% --------------------------------------------------------------------
    % Compute expectations taken wrt stationary distribution of baseline
    % ---------------------------------------------------------------------
    RA_with_BC.Empc =  HA_with_BC.pmf(:)' * RA_with_BC.mpcs(:);
    HA.Empc = HA_with_BC.pmf(:)' * HA.mpcs(:);

    %% --------------------------------------------------------------------
    % Decomposition
    % ---------------------------------------------------------------------
    for ia = 1:numel(p_baseline.abars)
        threshold = p_baseline.abars(ia);

        % Precompute P(b<threshold)
        p_htm = 0;
        for ii = 1:n_het
            p_htm = p_htm + HA_with_BC.cdf_interp{ii}(threshold);
        end

        % Term 1: RA MPC
        decomp.term1(ia) = RA.mpc;

        % Term 2: HtM effect
        mpc_integral = 0;
        for ii = 1:n_het
            mpc_integral = mpc_integral...
                + HA_with_BC.mpc_integral_interp{ii}(threshold);
        end
        decomp.term2(ia) = mpc_integral - RA.mpc * p_htm;

        % Term 3: Borrowing constraint effect
        mpc_integral = RA_with_BC.Empc;
        for ii = 1:n_het
            mpc_integral = mpc_integral...
                - RA_with_BC.mpc_integral_interp{ii}(threshold);
        end
        decomp.term3(ia) = mpc_integral - RA.mpc * (1-p_htm);

        % Term 4: Income risk effect
        mpc_integral = HA.Empc;
        for ii = 1:n_het
            mpc_integral = mpc_integral...
                - HA.mpc_integral_interp{ii}(threshold);
        end
        decomp.term4(ia) = mpc_integral - RA.mpc * (1-p_htm);

        % Term 5: Interaction
        decomp.term5(ia) = HA_with_BC.Empc - decomp.term1(ia)...
            - decomp.term2(ia) - decomp.term3(ia) - decomp.term4(ia);
    end
end

function initialized_decomp = initialize_to_nan(params)
	nan_vec = NaN(1, numel(params.abars));
	initialized_decomp.term1 = nan_vec;
	initialized_decomp.term2 = nan_vec;
	initialized_decomp.term3 = nan_vec;
	initialized_decomp.term4 = nan_vec;
	initialized_decomp.term5 = nan_vec;
end