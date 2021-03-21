function gini = direct_gini(level,distr)
	% Computes the gini index from a discrete probability
	% distribution.
	%
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

    % Sort distribution and levels by levels
    sorted = sortrows([level(:),distr(:)]);
    level_sort = sorted(:,1);
    dist_sort  = sorted(:,2);
    S = [0;cumsum(dist_sort .* level_sort)];
    gini = 1 - dist_sort' * (S(1:end-1)+S(2:end)) / S(end);
end