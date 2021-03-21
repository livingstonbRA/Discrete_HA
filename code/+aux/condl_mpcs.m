function mpcs_a = condl_mpcs(mpcs_states, pmf, pmf_a)
	% Computes mean MPCs over the asset space.

	na = size(mpcs_states, 1);
	mpcs_states = reshape(mpcs_states, na, []);
	pmf = reshape(pmf, na, []);

	mpcs_a = sum(mpcs_states .* pmf, 2)...
		./ pmf_a;

	% Use unweighted mean where marginal pmf is very small
	pmf_a_small = pmf_a < 1e-9;
	mpcs_a(pmf_a_small) = mean(mpcs_states(pmf_a_small,:), 2);
end