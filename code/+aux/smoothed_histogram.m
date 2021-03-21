function [bins, vals] = smoothed_histogram(agrid, pmf, nbins, amax)
	a_cdf = cumsum(pmf);
	[agrid_u, iu] = unique(agrid, 'last');
	a_cdf = a_cdf(iu);

	cdf_interp = griddedInterpolant(agrid_u, a_cdf, 'pchip', 'nearest');

	amin = agrid(1);
	spacing = (amax - amin) / nbins;
	bins = 1:nbins;
	bins = amin + bins * spacing;

	vals = zeros(nbins, 1);
	for ibin = 1:nbins
		bin_start = bins(ibin) - spacing;
		bin_end = bins(ibin);

		P_lt_start = cdf_interp(bin_start);

		if ibin < nbins
			P_lt_end = cdf_interp(bin_end);
		else
			P_lt_end = 1;
		end

		vals(ibin) = P_lt_end - P_lt_start;
	end

	bins = [amin, bins];
end