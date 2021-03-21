classdef Grid < handle
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties (SetAccess = private)
		% Grids for cash-on-hand, savings, and assets
		x;
		s;
		a;

		% Either 'EGP' or 'DST'
		gtype;

		% Cash-on-hand grid sizes
		nx;
		nx_neg;

		% Income grid sizes
		nyP;
		nyF;

		% Position of zero in the x-grid
		i0;

		% Returns matrix
		R_broadcast;
	end

	methods
		function obj = Grid(params, income, gtype)
			obj.gtype = gtype;

			if strcmp(gtype,'EGP')
                obj.nx = params.nx;
                obj.nx_neg = params.nx_neg;
			elseif strcmp(gtype,'DST')
				obj.nx = params.nx_DST;
				obj.nx_neg = params.nx_neg_DST;
			end

			if (params.borrow_lim<0) && (obj.nx_neg==0)
				error("Must have nx_neg > 0 since borrow_lim < 0")
			end

			obj.R_broadcast = reshape(params.R,...
				[1 1 1 numel(params.R)]);

			obj.create_sgrid(params);
			obj.create_xgrid(params, income);
			obj.create_norisk_xgrid(params, income);
            obj.create_agrid(params);

            obj.i0 = obj.nx_neg + 1;
		end

		function obj = create_sgrid(obj, params)
		    soft_constraint = 0;
			if obj.nx_neg > 0
				neg_midpt = (params.borrow_lim + soft_constraint) / 2;
				pts_bottom = round(obj.nx_neg/2);
				pts_top = obj.nx_neg - pts_bottom + 2;

				% Portion of savings grid closest to borrowing limit
				savgrid_neg_low = create_curved_grid(...
					params.borrow_lim, neg_midpt, pts_bottom,...
					params.xgrid_par_neg, false);

				% Portion of savings grid closest to soft constraint
				savgrid_neg_high = create_curved_grid(...
					neg_midpt, soft_constraint, pts_top,...
					params.xgrid_par_neg, true);

				% Combine while deleting repetitions
				savgrid_neg = [savgrid_neg_low; savgrid_neg_high(2:end-1)];
			else
				savgrid_neg = [];
			end

			nx_pos = obj.nx - obj.nx_neg;
            
            % Note to self: Why am I using the alt grid here??
            savgrid_pos = create_curved_grid_alt(soft_constraint,...
            	params.xmax, params.xgrid_term1wt,...
            	params.xgrid_term1curv, params.xgrid_par, nx_pos);

			savgrid = [savgrid_neg; savgrid_pos];

		    obj.s.vec = savgrid;
		    obj.s.matrix = repmat(savgrid, [1 params.nyP params.nyF]);
		end

		function obj = create_xgrid(obj, params, income)
			minyT = reshape(min(income.netymat, [], 2), [1 params.nyP params.nyF]);
			xgrid = obj.R_broadcast .* obj.s.matrix + minyT;
			if numel(params.r) == 1
				xgrid = repmat(xgrid, [1 1 1 params.nb]);
			end

			obj.x.matrix = xgrid;
		end

		function obj = create_norisk_xgrid(obj, params, income)
			xmatrix_norisk = obj.R_broadcast .* obj.s.vec + income.meannety1;
			xmatrix_norisk = reshape(xmatrix_norisk, obj.nx, []);

			if numel(params.r) == 1
				xmatrix_norisk = repmat(xmatrix_norisk, [1 params.nb]);
			end

			obj.x.matrix_norisk = xmatrix_norisk;
		end

		function obj = create_agrid(obj, params)
			% For returns heterogeneity:
			% agrid vector must extend as far downward as the lowest
			% possible assets, but allow for points close to soft constraint
			% as well
			R_selected = max(params.R) * (obj.s.vec < 0) ...
				+ min(params.R) * (obj.s.vec >= 0);
			obj.a.vec = R_selected .* obj.s.vec;

		 	obj.a.matrix = repmat(obj.a.vec, [1, params.nyP, params.nyF, params.nb]);
		end
	end

end

function cgrid = create_curved_grid(lowval, highval, npts,...
	curvature, reversed)
	cgrid = linspace(0, 1, npts)';
    cgrid = cgrid .^ (1/curvature);
    if reversed
    	cgrid = 1 - flip(cgrid);
    end
    cgrid = lowval + (highval - lowval) * cgrid;
end

function vgrid = create_curved_grid_alt(vmin, vmax,...
	term1_wt, term1_curv, curv, npts)
	vgrid = linspace(0, 1, npts)';
	vgrid = term1_wt * vgrid .^ (1 / term1_curv) ...
		+ vgrid .^ (1 / curv);
	vgrid = vgrid / vgrid(end);
	vgrid = vmin + (vmax - vmin) * vgrid;
end