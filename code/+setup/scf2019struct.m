function scf = scf2019struct()
	% Returns a structure containing aggregate statistics from the 2019 SCF.
	% HtM defined as household with b < y / 6
    % PHtM defined as household with (a + b) < y / 6

    scf = struct();
    scf.quarterly_earnings = 79181 / 4;
    scf.annual_earnings = 79181;
    scf.median_totw = 1.54;
    scf.median_liqw = 0.05;
    scf.htm = 0.39;
    scf.phtm = 0.135;
end