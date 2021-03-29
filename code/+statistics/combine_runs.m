
% This script combines .mat files named variablesX.mat into
% an excel spreadsheet

clear

[~, currdir] = fileparts(pwd());
if ~strcmp(currdir, 'Discrete_HA')
    msg = 'The user must cd into the Discrete_HA directory';
    bad_dir = MException('Discrete_HA:master', msg);
    throw(bad_dir);
end

% Check if code is running locally or on the server
taskid = str2num(getenv('SLURM_ARRAY_TASK_ID'));
running_on_server = ~isempty(taskid);

if ~running_on_server
    taskid = 1;
end

outdir = fullfile('output', sprintf('tables%d', taskid));
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

addpath('code');

% Read .mat files into a cell array
ind = 0;
for irun = 1:999
    fname = sprintf('variables%d.mat', irun);
    fpath = fullfile('output', fname);
    if exist(fpath,'file')
        ind = ind + 1;

        S = load(fpath);
        params(ind) = S.Sparams;
        results(ind) = S.results;
        stats{ind} = S.results.stats;
    end
end

if (ind == 0)
    error('No mat files found')
end

for ip = 1:ind
    if params(ip).freq == 1
        baseind = find(ismember({params.name}, {'Annual'}));
    else
        baseind = find(ismember({params.name}, {'Quarterly'}));
    end

    if isempty(baseind)
        baseind = ip;
    end

    p0 = params(baseind);
    p1 = params(ip);
    stats0 = stats{baseind};
    stats1 = stats{ip};
    cdecomp = statistics.ComparisonDecomp(p0, p1, stats0, stats1);

    mpcs0 = reshape(stats0.mpcs(5).mpcs_1_t{1},...
        p0.nx_DST, []);
    mpcs1 = reshape(stats1.mpcs(5).mpcs_1_t{1},...
        p1.nx_DST, []);
    cdecomp.perform_decompositions(mpcs0, mpcs1);
    decomps_baseline(ip) = cdecomp.results;
end

ctimename = sprintf('continuous_time_baseline%d.mat', taskid);
ctimepath = fullfile('input', ctimename);
ctimeresults = tables.read_continuous_time_results(ctimepath);

tables.TexTables.save_baselines_tables(params, results, outdir, 'ctimeresults', ctimeresults);

for ip = 3:11
	tables.TexTables.save_experiment_table(params, results, decomps_baseline, outdir, ip);
end