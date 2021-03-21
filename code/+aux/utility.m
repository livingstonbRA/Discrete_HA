function u = utility(risk_aver,c)
	if isequal(size(risk_aver),size(c))
		u = zeros(size(risk_aver));
		Ilog_u = (risk_aver == 1);

		u(Ilog_u) = log(c(Ilog_u));
		u(~Ilog_u) = (c(~Ilog_u) .^(1-risk_aver(~Ilog_u) )-1)./(1-risk_aver(~Ilog_u) );
	elseif numel(risk_aver) == 1
		if risk_aver == 1
			u = log(c);
		else
			u = (c .^(1-risk_aver) - 1) ./ (1 - risk_aver);
		end
	end
end