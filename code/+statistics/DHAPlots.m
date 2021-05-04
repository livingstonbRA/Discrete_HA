classdef DHAPlots
    properties (Constant)
        fontsize = 13;
    end

    methods (Static)
        function consumption_wealth_overlay(stats, p, varargin)
            parser = inputParser;
            addOptional(parser, 'x_lb', 0);
            addOptional(parser, 'x_ub', 10);
            addOptional(parser, 'nbins', 40);
            addOptional(parser, 'curve_variable', 'yP');
            addOptional(parser, 'curve_indices', []);
            parse(parser, varargin{:});

            options = parser.Results;
            
            xvars = stats.xgrid_variables;
            
            % Cash-on-hand histogram
            sorted_mat = sortrows([xvars.xvals(:), xvars.xdist(:)]);
            xvals = sorted_mat(:,1);
            xdist = sorted_mat(:,2);

            [xvals, ixunique] = unique(xvals);
            xid = zeros(numel(xvars.xvals), 1);
            xid(ixunique) = 1;
            xid = cumsum(xid);
            xdist = accumarray(xid, xdist);
            
            [edges, counts] = smoothed_histogram(xvals, xdist, options.nbins, options.x_ub);
            
            % Consumption policy functions
            if strcmp(options.curve_variable, 'yP')
                if isempty(options.curve_indices)
                    iy_median = ceil(p.nyP / 2);
                    options.curve_indices = [1, iy_median, p.nyP];
                    con_x = cell(1, 3);
                    for ifig = 1:3
                        iyP = options.curve_indices(ifig);
                        con_x{ifig} = sortrows([xvars.xvals(:,iyP), xvars.con_x(:,iyP)]);
                    end
                    labels = {'$y_{high}$', '$y_{mid}$', '$y_{low}$'};
                    figure_order = [3, 2, 1];
                end
            elseif ismember({options.curve_variable}, {'beta', 'crra'})
                iyP = ceil(p.nyP / 2);
                if isempty(options.curve_indices)
                    iz_median = ceil(p.nz / 2);
                    options.curve_indices = [1, iz_median, p.nz];
                    
                end

                nfigs = numel(options.curve_indices);
                con_x = cell(1, nfigs);

                for ifig = 1:nfigs
                    iz = options.curve_indices(ifig);
                    xvals_iz = xvars.xvals(:,iyP,1,iz);
                    cvals_z = xvars.con_x(:,iyP,1,iz);
                    con_x{ifig} = sortrows([xvals_iz(:), cvals_z(:)]);
                end
            
                if strcmp(options.curve_variable, 'beta')
                    labels = {'$\beta_{high}$', '$\beta_{mid}$', '$\beta_{low}$'};
                    figure_order = [3, 2, 1];
                elseif strcmp(options.curve_variable, 'crra')
                    labels = {'$crra_{high}$', '$crra_{mid}$', '$crra_{low}$'};
                    figure_order = [3, 2, 1];
                end
            end
            
            % Plots
            close all
            ax = axes();

            hold(ax, 'on')
            for ifig = figure_order
                plot(con_x{ifig}(:,1), con_x{ifig}(:,2), 'Parent', ax)
            end
            hold(ax, 'off')

            ylabel('$c_t$', 'FontWeight', 'bold', 'Interpreter', 'latex',...
                'Rotation', 0)

            yyaxis right
            wealth_hist = histogram('Parent', ax, 'BinEdges', edges, 'BinCounts', counts);
            ylabel('$f_x$', 'FontWeight', 'bold', 'Interpreter', 'latex',...
                'Rotation', 0)
            
            xlabel('$x_t$', 'FontWeight', 'bold', 'Interpreter', 'latex')
            legend(labels, 'Interpreter','latex', 'Location', 'southeast')
            legend('boxoff')

            format_plot(ax);
            xlim([0, options.x_ub])
            
            yyaxis left
            ylab = get(ax, 'ylabel');
            ylab.Position

            if strcmp(options.curve_variable, 'beta')
                ylab.Position(1) = -1;
            else
                ylab.Position(1) = -1;
            end
        end

        function mpc_function(stats, varargin)
            parser = inputParser;
            addOptional(parser, 'mpc_type', 'period');
            parse(parser, varargin{:});

            options = parser.Results;

            if strcmp(options.mpc_type, 'period')
                mpcs = stats.mpcs(5).avg_s_t(1,1:9);
                label = '$MPC_t$';
            else
                mpcs = cumsum(stats.mpcs(5).avg_s_t(1,1:9));
                label = '$MPC_t$';
            end 

            close all
            ax = axes();
            plot(0:8,mpcs);
            xlabel('$t$', 'FontWeight', 'bold', 'Interpreter', 'latex')
            ylabel(label, 'FontWeight', 'bold', 'Interpreter', 'latex',...
                'Rotation', 0)

            ylab = get(ax, 'ylabel');
            ylab.Position(1) = -1.5;

            pos=get(gca,'position');
            xshift = 0.08;
            pos(1) = pos(1) + xshift;
            % pos(3) = pos(3) - xshift / 2;
            set(gca,'position',pos);  % write the new values

            ax = format_plot(ax);
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

function ax = format_plot(ax)
    set(gcf,'color','w');
    set(ax, 'FontSize', statistics.DHAPlots.fontsize);
end