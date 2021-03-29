classdef Calibrator < handle
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	properties
		% Calibrator options
		options;

		% Names of parameters to vary (in Params instance), cell array
		variables;

		% Names of target variables (in Statistics instance), cell array
		target_names;

		% Values to match for target variables, numeric vector
		target_values;

		% Bounds on parameters
		lbounds = [];
		ubounds = [];

		% Scale factors for target variable devations, optional
		fscale = [];

		% Number of parameters/target variables
		n;

		% Initial parameter values
		x0;

		% Computation
		iter = 1;

		% Function handle that calls code to solve the model
		main_handle;

		% Function handle that accesses equilibrium values of target variables
		results_handle;

		% Function handle to pass to nonlinear solver
		solver_handle;
	end

	methods
		function obj = Calibrator(params, variables,...
			target_names, target_values)
			obj.construct_options_struct(params);

			obj.variables = variables;
			obj.n = numel(variables);

			obj.target_names = target_names;
			obj.target_values = target_values;

			obj.fscale = ones(obj.n, 1);

			x0_1 = zeros(1, obj.n);
			for i_var = 1:obj.n
				x0_1(i_var) = params.(obj.variables{i_var});
			end
			
			obj.x0 = {x0_1};
		end

		function construct_options_struct(obj, params)
			obj.options = struct();
		end

		function set_fscale(obj, fscale_in)
			obj.fscale = fscale_in(:);
		end

		function set_param_bounds(obj, varargin)
			nv = numel(varargin);
			assert(nv==obj.n, "Too many or too few bounds provided");

			for ii = 1:nv
				var_bounds = varargin{ii};
				obj.lbounds(ii) = var_bounds(1);
				obj.ubounds(ii) = var_bounds(2);
			end

			for ix0 = 1:numel(obj.x0)
				x0 = obj.x0{ix0};
				for ii = 1:obj.n
					if x0(ii) < obj.lbounds(ii)
						x0(ii) = (5*obj.lbounds(ii) + obj.ubounds(ii))/6;
					elseif x0(ii) > obj.ubounds(ii)
						x0(ii) = (obj.lbounds(ii) + 5*obj.ubounds(ii))/6;
                    end
                    obj.x0{ix0} = x0;
				end
			end
		end

		function set_handle(obj, p)
			if isempty(obj.main_handle)
				error('Derived class has not initialized main_handle')
			end

			obj.solver_handle = @(x) obj.fn_handle(x, p);
		end

		function dv = fn_handle(obj, x, current_params)
			obj.turn_off_param_options(current_params);
			quiet = true;

			for i_var = 1:obj.n
				current_params.set(obj.variables{i_var}, x(i_var), quiet);
			end
			results = obj.main_handle(current_params);

			fprintf('  -- function evaluation %d --\n    evaluated at: ', obj.iter)
			v = zeros(obj.n, 1);
			for i_var = 1:obj.n
				v(i_var) = obj.get_results_value(results, obj.target_names{i_var});
				fprintf('%s = %g', obj.variables{i_var}, x(i_var))
				if i_var < obj.n
					fprintf(", ")
				end
			end
			dv = v(:) - obj.target_values(:);
			dv = dv .* obj.fscale;

			fprintf('\n    target variables: ')
			for i_var = 1:obj.n
				fprintf('%s = %g', obj.target_names{i_var}, v(i_var))
				if i_var < obj.n
					fprintf(", ")
				end
			end
			fprintf('\n    norm: %f\n', norm(dv))

			obj.reset_param_options(current_params);
			obj.iter = obj.iter + 1;
		end

		function value = get_results_value(obj, results, variable_name)
			value = nan;
		end

		function solver_args = get_args(obj)
			x0 = obj.x0{1};
			if ~isempty(obj.ubounds) && ~isempty(x0)
				solver_args = {x0, obj.lbounds, obj.ubounds};
			elseif ~isempty(x0)
				solver_args = {x0};
			else
				solver_args = {};
			end
		end

		function turn_off_param_options(obj, params)
			quiet = true;
			props = fields(obj.options);
			for ip = 1:numel(props)
				params.set(props{ip}, 0, quiet);
			end
		end

		function reset_param_options(obj, params)
			quiet = true;
			props = fields(obj.options);
			for ip = 1:numel(props)
				params.set(props{ip}, obj.options.(props{ip}), quiet);
			end
		end
	end
end
