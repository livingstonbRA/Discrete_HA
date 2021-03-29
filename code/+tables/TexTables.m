classdef TexTables
    
    properties (Constant)
        table_includes = {
            {'Baseline'}
            {'Baseline'}
            {'Q1a'}
            {'Q1b'}
            {'Q2'}
            {'Q3'}
            {'Q4'}
            {'Q5'}
            {'Q6'}
            {'Q7'}
            {'Q8'}
        };
    end

    methods (Static)
        function save_baselines_tables(params_in, results, dirpath, varargin)
            for tableno = [1, 2]
	            for panelname = {'header', 'A', 'B', 'C', 'D'}
	            	if (tableno == 1)
	            		panelobj = tables.TexTables.table1panel(params_in, results, panelname{:}, varargin{:});
	            	elseif (tableno == 2)
	            		panelobj = tables.TexTables.table2panel(params_in, results, panelname{:}, varargin{:});
	            	end

	            	if strcmp(panelname{:}, 'header')
	            		panelfname = sprintf('table%d_header.xlsx', tableno);
	            	else
	            		panelfname = sprintf('table%d_panel%s.xlsx', tableno, panelname{:});
	            	end

	            	panelfpath = fullfile(dirpath, panelfname);
	            	writetable(panelobj, panelfpath, 'WriteRowNames', true);
	            end
	        end
        end

        function save_experiment_table(params_in, results, comparison_decomps, dirpath, tableno)
            header = tables.TexTables.experiment_table_header(params_in, results, tableno);
            headerpath = fullfile(dirpath, sprintf('table%d_header.xlsx', tableno));
            writetable(header, headerpath, 'WriteRowNames', true);

            for panelname = {'A', 'A2', 'B', 'C', 'D'}
            	if ismember(panelname{:}, {'A', 'A2'})
            		panelobj = tables.TexTables.experiment_table_panel(...
            			params_in, comparison_decomps, panelname{:}, tableno);
            	else
            		panelobj = tables.TexTables.experiment_table_panel(...
            			params_in, results, panelname{:}, tableno);
            	end
            	panelfname = sprintf('table%d_panel%s.xlsx', tableno, panelname{:});
            	panelfpath = fullfile(dirpath, panelfname);
            	writetable(panelobj, panelfpath, 'WriteRowNames', true);
            end
        end

        function table_out = table1panel(params_in, results, panel, varargin)
            switch panel
	            case 'header'
	            	get_stats = @(x) {
	            		x.stats.mpcs(5).quarterly
	                    x.stats.mpcs(5).annual
	                    x.stats.beta_A
	                  };
	            case 'A'
		            get_stats = @(x) {...
		            	x.stats.mean_gross_y_annual
                        x.stats.std_log_gross_y_annual
                        x.stats.std_log_net_y_annual
	                };
	            case 'B'
	            	get_stats = @(x) {...
	            		x.stats.mean_a
                        x.stats.sav0
                        x.stats.constrained{1}
                        x.stats.constrained_dollars{1}
                        x.stats.constrained_dollars{2}
                        x.stats.constrained_dollars{3}
                        x.stats.constrained_dollars{4}
                        x.stats.a_lt_ysixth
                        x.stats.a_lt_ytwelfth
                        x.stats.wpercentiles{1}
                        x.stats.wpercentiles{2}
                        x.stats.wpercentiles{3}
                        x.stats.wpercentiles{5}
                        x.stats.wpercentiles{7}
                        x.stats.wpercentiles{8}
                        x.stats.w_top10share
                        x.stats.w_top1share
                        x.stats.wgini
                    };
                case 'C'
                	get_stats = @(x) {...
                		x.stats.mpcs(4).annual
                        x.stats.mpcs(6).annual
                        x.stats.mpcs(4).quarterly
                        x.stats.mpcs(6).quarterly
                    };
                case 'D'
                	get_stats = @(x) {...
                		x.stats.mpcs(1).annual
                        x.stats.mpcs(2).annual
                        x.stats.mpcs(3).annual
                        x.stats.mpcs(1).quarterly
                        x.stats.mpcs(2).quarterly
                        x.stats.mpcs(3).quarterly
                    };
                otherwise
                	error("Invalid panel entry")
            end

            [stats_array, names_array] = stack_results(1, get_stats, params_in, results, varargin{:});
            table_out = make_table(stats_array, names_array);
        end

        function table_out = table2panel(params_in, results, panel, varargin)
            decomp_norisk_get_stats_fn = @(x, k)  {
            	x.stats.decomp_norisk.term2(k)
                x.stats.decomp_norisk.term3(k)
                x.stats.decomp_norisk.term4(k)
            };

            switch panel
	            case 'header'
	            	get_stats = @(x) {
	            		x.stats.mpcs(5).oneperiod
	            		x.stats.decomp_norisk.term1_pct
	            	};
	            case 'A'
	                get_stats = @(x) decomp_norisk_get_stats_fn(x, 1);
	            case 'B'
	            	get_stats = @(x) decomp_norisk_get_stats_fn(x, 2);
	           	case 'C'
	           		get_stats = @(x) decomp_norisk_get_stats_fn(x, 3);
	            case 'D'
	            	get_stats = @(x) {
	            		x.stats.decomp_RA.Em1_less_mRA
	                    x.stats.decomp_RA.term1
	                    x.stats.decomp_RA.term2
	                    x.stats.decomp_RA.term3
	                };
	           	otherwise
	           		error("Invalid panel selection")
           	end

            [stats_array, names_array] = stack_results(2, get_stats, params_in, results, varargin{:});
            table_out = make_table(stats_array, names_array);
        end

        function table_out = experiment_table_header(params_in, results, tableno)
            indices = filter_param_group(params_in, tables.TexTables.table_includes{tableno});

            import statistics.Statistics.sfill

            for ii = 1:numel(indices)
                ip = indices(ii);
                if (tableno == 10) || (tableno == 11)
                    if ~isempty(params_in(ip).tex_header_values)
                        tex_vals = params_in(ip).tex_header_values{1};
                        variable_values = {
                            sfill(string(tex_vals.description), 'Description')
                        };
                    else
                        variable_values = {
                            sfill(nan, 'Description')
                        };
                    end
                elseif (tableno == 9)
                    if ~isempty(params_in(ip).tex_header_values)
                        tex_vals = params_in(ip).tex_header_values{1};
                        variable_values = {
                            sfill(tex_vals.r, 'r')
                        };
                    else
                        variable_values = {
                            sfill(nan, 'r', 3)
                        };
                    end
                elseif (tableno == 8)
                    if ~isempty(params_in(ip).tex_header_values)
                        tex_vals = params_in(ip).tex_header_values{1};
                        variable_values = {
                            sfill(tex_vals.tempt, 'Temptation')
                        };
                    else
                        variable_values = {
                            sfill(nan, 'Temptation', 3)
                        };
                    end
                elseif (tableno == 6) || (tableno == 7)
                    if ~isempty(params_in(ip).tex_header_values)
                        tex_vals = params_in(ip).tex_header_values{1};
                        variable_values = {
                            sfill(string(tex_vals.riskaver), 'Risk aversion')
                            sfill(tex_vals.ies, 'IES')
                        };
                    else
                        variable_values = {
                            sfill(nan, 'Risk aversion', 2)
                            sfill(nan, 'IES', 3)
                        };
                    end
                elseif tableno == 5
                    if ~isempty(params_in(ip).tex_header_values)
                        tex_vals = params_in(ip).tex_header_values{1};
                        variable_values = {
                            sfill(tex_vals.pswitch, 'Switch probability', 2)
                            sfill(tex_vals.width, 'Spacing', 3)
                        };
                    else
                        variable_values = {
                            sfill(nan, 'Switch probability', 2)
                            sfill(nan, 'Spacing', 3)
                        };
                    end
                elseif tableno == 4
                    variable_values = {};
                elseif tableno == 3
                    if ~isempty(params_in(ip).tex_header_values)
                        tex_vals = params_in(ip).tex_header_values{1};
                        variable_values = {sfill(tex_vals.value, 'Value', 2)};
                    else
                        variable_values = {sfill(nan, 'Value', 2)};
                    end
                end
                if tableno == 8
                    statistics{ii} = {  results(ip).stats.mpcs(5).annual
                                        results(ip).stats.beta_A
                                      };
                else
                    statistics{ii} = {  results(ip).stats.mpcs(5).quarterly
                                        results(ip).stats.mpcs(5).annual
                                        results(ip).stats.beta_A
                                      };
                end
                statistics{ii} = [variable_values(:); statistics{ii}];
                names{ii} = params_in(ip).tex_header;
            end

            table_out = make_table(statistics, names, 'experiment', true);
        end

        function table_out = experiment_table_panel(params_in, variables, panel, tableno)
            switch panel
	            case 'A'
	            	get_stats = @(x) {
	            		x.Em1_less_Em0
                        x.term1
                        x.term2
                        x.term2a(2)
                        x.term2b(2)
                        x.term3
	            	};
	           	case 'A2'
	           		get_stats = @(x) {
	           			x.term1_pct
                        x.term2_pct
                        x.term2a_pct(2)
                        x.term2b_pct(2)
                        x.term3_pct
	           		};
	           	case 'B'
	           		get_stats = @(x) {
	           			x.stats.mean_a
	                    x.stats.sav0
	                    x.stats.constrained{1}
	                    x.stats.constrained_dollars{1}
	                    x.stats.constrained_dollars{2}
	                    x.stats.constrained_dollars{3}
	                    x.stats.constrained_dollars{4}
	                    x.stats.a_lt_ysixth
	                    x.stats.a_lt_ytwelfth
	                    x.stats.w_top10share
	                    x.stats.w_top1share
	                    x.stats.wgini
	           		};
	           	case 'C'
                    if tableno == 8
                        get_stats = @(x) {
                            x.stats.mpcs(4).annual
                            x.stats.mpcs(6).annual
                        };
                    else
    	           		get_stats = @(x) {
    	           			x.stats.mpcs(4).quarterly
                            x.stats.mpcs(6).quarterly
                        };
                    end
                case 'D'
                    if tableno == 8
                        get_stats = @(x) {
                            x.stats.mpcs(1).annual
                            x.stats.mpcs(2).annual
                            x.stats.mpcs(3).annual
                        };
                    else
    	           		get_stats = @(x) {
    	           			x.stats.mpcs(1).quarterly
                            x.stats.mpcs(2).quarterly
                            x.stats.mpcs(3).quarterly
                        };
                    end
            end

            [stats_array, names_array] = stack_results(tableno, get_stats, params_in, variables, 'experiment', true);
            table_out = make_table(stats_array, names_array, 'experiment', true);
        end
    end 
end

function indices = filter_param_group(params_in, includes)
    indices = [];
    jj = 1;
    for ii = 1:numel(params_in)
        keep = false;
        for kk = 1:numel(params_in(ii).group)
            if ismember(params_in(ii).group{kk}, includes)
                keep = true;
                break;
            end
        end

        if keep
            indices(jj) = ii;
            jj = jj + 1;
        end
    end
end

function [stats_array, names_array] = stack_results(tableno, fn_handle, params_in, main_results, varargin)
    parser = inputParser;
    addOptional(parser, 'experiment', false);
    addOptional(parser, 'ctimeresults', []);
    parse(parser, varargin{:});
    experiment = parser.Results.experiment;
    ctimeresults = parser.Results.ctimeresults;

    indices = filter_param_group(params_in, tables.TexTables.table_includes{tableno});
    n = numel(indices);
    stats_array = cell(n, 1);
    names_array = cell(n, 1);
    for ii = 1:n
        ip = indices(ii);
        stats_array{ii} = fn_handle(main_results(ip));

        if experiment & ~isempty(params_in(ip).tex_header)
            names_array{ii} = params_in(ip).tex_header;
        else
            names_array{ii} = params_in(ip).name;
        end
    end

    if ~isempty(ctimeresults)
        stats_array{n+1} = fn_handle(ctimeresults);
        names_array{n+1} = 'Continuous Time';
    end
end

function table_out = make_table(statistics, names, varargin)
    parser = inputParser;
    addOptional(parser, 'experiment', false);
    parse(parser, varargin{:});
    experiment = parser.Results.experiment;

    n = numel(names);
    for ii = 1:n
        vars{ii} = get_values(statistics{ii});

        if experiment
            varnames{ii} = char(names{ii} + string(sprintf('__v%d__', ii)));
        else
            varnames{ii} = names{ii};
        end
    end

    vars{n+1} = get_precision(statistics{1});
    varnames{n+1} = 'decimals';
    rownames = get_names(statistics{1});
    table_out = table(vars{:}, 'RowNames', rownames(:), 'VariableNames', varnames(:));
end

function values = get_values(entries)
    values = {};
    for ii = 1:numel(entries)
        values{ii} = entries{ii}.value;
    end
    values = values(:);
end

function decimals = get_precision(entries)
    decimals = {};
    for ii = 1:numel(entries)
        decimals{ii} = entries{ii}.decimals;
    end
    decimals = decimals(:);
end

function names = get_names(entries)
    names = {};
    for ii = 1:numel(entries)
        names{ii} = entries{ii}.tex_label;
    end
    names = names(:);
end