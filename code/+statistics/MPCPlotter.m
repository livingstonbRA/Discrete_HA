classdef MPCPlotter < handle

	properties
		p;
		grids;
		dims;
		mpcs;

		lims = [0, 10, -inf, inf];
		fig;
		ax;
		tile;

		mpcs_plotted = false;

		fontsize = 12;
		show_grid = 'off';
	end

	methods
		function obj = MPCPlotter(params, grids, mpcs)
			obj.p = params;
			obj.dims = [size(grids,1), obj.p.nyP, obj.p.nyF, obj.p.nz];

			obj.grids = grids;
			obj.mpcs = mpcs;
		end

		function [ax_main, ax_window] = create_mpcs_plot_yPs(...
			obj, yP_indices, zoomed_window, shock_size)

			obj.mpcs = reshape(obj.mpcs, obj.dims);

			obj.fig = figure();
			ax_main = axes(obj.fig);
			legend('hide');
			xlabel("Wealth (ratio to mean annual income)");
			ylabel(sprintf("MPC out of %g", shock_size));
			ax_main = obj.apply_formatting(ax_main);

			if zoomed_window
				ax_window = add_window(obj.fig, ax_main);
				ax_window = obj.apply_formatting(ax_window);
			end


			for ii = 1:numel(yP_indices)
			    iyP = yP_indices(ii);
			    mpcs_yP = obj.mpcs(:,iyP,1,1);
			    agrid = obj.grids(:,iyP,1,1);

			    obj.plot_mpc_function(ax_main, mpcs_yP, agrid);
			    if zoomed_window
			        obj.plot_mpc_function(ax_window, mpcs_yP, agrid);
			    end
			end

			yP_labels = {'Low yP', 'High yP'};
			legend(ax_main, yP_labels);
			legend(ax_main, 'show');
			box(ax_main.Legend, 'off');
		end

		function [ax_main, ax_window] = create_mpc_plot_shocks(...
			obj, yP_index, zoomed_window, ishocks)

			obj.fig = figure();
			ax_main = axes(obj.fig);
			legend('hide');
			xlabel("Wealth (ratio to mean annual income)");
			ylabel("MPC");
			ax_main = obj.apply_formatting(ax_main);

			if zoomed_window
				ax_window = add_window(obj.fig, ax_main);
				ax_window = obj.apply_formatting(ax_window);
			end

			shock_labels = {};
			for ii = 1:numel(ishocks)
				ishock = ishocks(ii);
			    mpcs_ishock = obj.mpcs{ii}(:,yP_index,1,1);
			    agrid = obj.grids(:,yP_index,1,1);

			    obj.plot_mpc_function(ax_main, mpcs_ishock, agrid);
			    if zoomed_window
			        obj.plot_mpc_function(ax_window, mpcs_ishock, agrid);
			    end

			    shock_labels{ii} = sprintf('Shock of %g', obj.p.shocks(ishock));
			end

			legend(ax_main, shock_labels);
			legend(ax_main, 'show');
			box(ax_main.Legend, 'off');
		end

		function plot_mpc_function(obj, parent_obj, mpcs, agrid)
			hold(parent_obj, 'on');
			plot(parent_obj, agrid, mpcs);
			hold(parent_obj, 'off');
		end

		function ax = add_median_wealth(obj, ax, median_wealth)
			hold(ax, 'on')
			wline = [0, 1];
			wvals = [median_wealth, median_wealth];
			plot(ax, wvals, wline, '--b');
		    ax.Legend.String{end} = 'Median wealth';
		    hold(ax, 'off')
		end

		function ax = apply_formatting(obj, ax)
			grid(ax, obj.show_grid);
			set(ax, 'FontSize', obj.fontsize);
		end
	end
end

function ax = add_window(parent_obj, main_ax)
	main_ax_corner = main_ax.Position(1:2);
	main_ax_width = main_ax.Position(3);
	main_ax_height = main_ax.Position(4);

	window_corner(1) = main_ax_corner(1) + 0.25 * main_ax_width;
	window_corner(2) = main_ax_corner(2) + 0.4 * main_ax_height;
	window_width = 0.5 * main_ax_width;
	window_height = 0.3 * main_ax_height;

	window_position = [window_corner, window_width, window_height];
	ax = axes('Position', window_position);
	box(ax, 'on');
end
