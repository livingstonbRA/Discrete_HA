classdef WealthPlotter

	properties
		p;
		dims;
		agrid;
		pmf_a;

		fig;

		fontsize = 12;
	end

	methods
		function obj = WealthPlotter(params, agrid, pmf_a)
			obj.p = params;
			obj.dims = [obj.p.nx_DST, obj.p.nyP, obj.p.nyF, obj.p.nb];

			obj.agrid = agrid;
			obj.pmf_a = pmf_a;
		end

		function [ax, wealth_hist] = create_histogram(obj, nbins, amax)
			if nargin < 3
				amax = max(obj.agrid);
			end

			obj.fig = figure();
			[edges, counts] = smoothed_histogram(obj.agrid, obj.pmf_a, nbins, amax);
			wealth_hist = histogram('Parent', obj.fig, 'BinEdges', edges, 'BinCounts', counts);

			bar_color = [0,0,0] + 0.5;
			wealth_hist.FaceColor = bar_color;
			wealth_hist.EdgeColor = bar_color;

			ax = gca;
			xlabel("Wealth (ratio to mean annual income)")
			ylabel("Probability density")
			set(ax, 'FontSize', obj.fontsize);
        end
	end
end

function [bins, vals] = smoothed_histogram(agrid, pmf, nbins, amax)
	a_cdf = cumsum(pmf);
	cdf_interp = griddedInterpolant(agrid, a_cdf,...
		'pchip', 'nearest');

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