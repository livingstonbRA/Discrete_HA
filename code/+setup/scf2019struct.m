function scf = scf2019struct()
	% Returns a structure containing aggregate statistics from the 2019 SCF.
	% HtM defined as household with b < y / 6
    % PHtM defined as household with (a + b) < y / 6

    scf = struct();
    scf.quarterly_earnings = 67131.733 / 4;
    scf.annual_earnings = 67131.733;
    scf.mean_totw = 275664.74; % 4.1
    scf.mean_liqw = 37709.078;
    scf.median_totw = 103380; % 1.54
    scf.median_liqw = 3100; % 0.046
    scf.htm = 0.4;
    scf.phtm = 0.1505;
end