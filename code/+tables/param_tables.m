function ptable =  param_tables(params, heterogeneity)

	param_names = {};

	param_names(1,:) = {'numeraire_in_dollars', 'Mean annual income ($)'};
	param_names(2,:) = {'nx', 'Num points, cash-on-hand grid'};
	param_names(3,:) = {'xmax', 'Max, cash-on-hand grid'};
	param_names(4,:) = {'xgrid_par', 'Curv, cash-on-hand grid (1 is linear)'};
	param_names(5,:) = {'borrow_lim', 'Borrowing limit'};
	param_names(6,:) = {'r', 'r'};
	param_names(7,:) = {'R', 'R'};
	param_names(8,:) = {'annuities', 'Annuity markets, off or on (0 or 1)'};
	param_names(9,:) = {'dieprob', 'One-period probability of death'};
	param_names(10,:) = {'invies', 'Inverse IES'};
	param_names(11,:) = {'risk_aver', 'Coeff rel risk aversion'};
	param_names(12,:) = {'beta0', 'One-period time discount factor (beta)'};
	param_names(13,:) = {'Bequests', 'Bequests, off or on (0 or 1)'};
	param_names(14,:) = {'nyT', 'Num points, transitory income'};
	param_names(15,:) = {'nyP', 'Num points, persistent income'};
	param_names(16,:) = {'nyF', 'Num points, fixed income'};
	param_names(17,:) = {'labtaxlow', 'Income tax rate (proportional)'};
	param_names(18,:) = {'labtaxhigh', 'Income tax rate (above threshold)'};
	param_names(19,:) = {'labtaxthreshpc', 'High income tax threshold (percentile)'};
	param_names(20,:) = {'savtax', 'Tax rate on savings'};
	param_names(21,:) = {'lumptransfer', 'One-period value of lump-sum transfer'};
	param_names(22,:) = {'nbeta', 'Num points, discount factor grid'};
	param_names(23,:) = {'betawidth', 'Distance between discount factor values'};
    param_names(24,:) = {'betagrid', 'One-period time discount factor(s)'};
	param_names(25,:) = {'nz', 'Num points, extra dim of heterogeneity (z-dim)'};
	param_names(26,:) = {'temptation', 'Temptation'};
    
	np = numel(params);
	nv = size(param_names, 1);
	values = cell(nv+1, np);
	variable_names = cell(1, np);
	row_names = cell(1, nv+1);

    row_names{1} = 'Description';
	for iv = 1:nv
		row_names{iv+1} = param_names{iv,2};
	end

	for ip = 1:np
        variable_names{ip} = sprintf('Run%d', ip);

        values{1,ip} = params(ip).name;
		for iv = 1:nv
			pname = param_names{iv,1};
            if strcmp(pname, 'betagrid')
                values{iv+1,ip} = num2str(heterogeneity(ip).(pname));
            elseif isnumeric(params(ip).(pname))
                values{iv+1,ip} = num2str(params(ip).(pname));
            else
                values{iv+1,ip} = params(ip).(pname);
            end
		end
    end
    
	ptable = cell2table(values, 'VariableNames', variable_names,...
		'RowNames', row_names);
end
