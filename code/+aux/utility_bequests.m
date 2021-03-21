function u = utility_bequests(bequest_curv, bequest_weight,...
	bequest_luxury,a)

	if bequest_curv == 1
        u = bequest_weight .* log(a + bequest_luxury);
    else
        u = bequest_weight .* ((a + bequest_luxury) ...
        	.^ (1-bequest_curv)-1)./(1-bequest_curv);
    end
end