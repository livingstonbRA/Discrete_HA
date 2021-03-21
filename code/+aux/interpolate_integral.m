function interpolant = interpolate_integral(...
	gridValues, integrandValues, pmf, is_sorted)
	% Returns an interpolant that interpolates to find the value of
	% int_0^{epsilon} values(a)g(a)da for a given epsilon.
	%
	% Evaluate the above integral by calling interpolant(epsilon).
	%
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	if nargin < 4
		is_sorted = false;
	end

	if ~is_sorted
		sortedInputs = sortrows(...
			[gridValues(:) integrandValues(:) pmf(:)]);
		gridSorted = sortedInputs(:,1);
		integrandSorted = sortedInputs(:,2);
		pmfSorted = sortedInputs(:,3);
	else
		gridSorted = gridValues;
		integrandSorted = integrandValues;
		pmfSorted = pmf;
	end

	integralValues= cumsum(integrandSorted .* pmfSorted);
	[gridUnique, uniqueInds] = unique(gridSorted, 'last');
	integralUnique = integralValues(uniqueInds);

	interpolant = griddedInterpolant(gridUnique, integralUnique,...
		'pchip', 'nearest');
end