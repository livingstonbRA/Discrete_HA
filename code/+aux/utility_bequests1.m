function u1 = utility_bequests1(bequest_curv, bequest_weight,...
	bequest_luxury, a)

	u1 = bequest_weight .* (a + bequest_luxury) .^ (-bequest_curv);
end