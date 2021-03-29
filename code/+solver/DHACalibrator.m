classdef DHACalibrator < solver.Calibrator
	% Brian Livingston, 2020
	% livingstonb@uchicago.edu

	methods
		function obj = DHACalibrator(params, variables,...
			target_names, target_values)
			% ARGUMENTS:
			%
			% params
			%	A Params object containing model parameters
			%
			% variables
			%	A cell array containing string vectors identifying the parameters
			%	to be calibrated, e.g. {'beta0'}.
			%
			% target_names
			%	A cell array containing the field names of the statistics to target.
			%	Each element should correspond to a field name of the results.direct
			%	structure. E.g {'median_a'}.
			%
			% target_values
			%	A numerical vector containing the target values for the targeted statistics.
			%	E.g. if the intent is to match median assets to 1.4, then target_names
			%	should include 'median_a' and target_values should include 1.4.

			obj = obj@solver.Calibrator(params, variables, target_names, target_values);
			
			obj.main_handle = @(curr_params) main(curr_params);
		end

		function construct_options_struct(obj, params)
			obj.options = struct();
			obj.options.MPCs = params.MPCs;
			obj.options.MPCs_news = params.MPCs_news;
			obj.options.MPCs_loan_and_loss = params.MPCs_loan_and_loss;
			obj.options.Simulate = params.Simulate;
			obj.options.DeterministicMPCs = params.DeterministicMPCs;
			obj.options.MakePlots = params.MakePlots;
			obj.options.SaveOutput = params.SaveOutput;
		end

		function value = get_results_value(obj, results, variable_name)
			value = results.stats.(variable_name).value;
		end
	end
end