
clear

cd D:\
cd GitHub/Discrete_HA/output
cd tables_7_26_21

tempfile newtable

local i = 51

import excel using "table`i'.xlsx", allstring
rename B run`i'
gen nid = _n
save newtable, replace


while `i' > 1 {
	local --i
	import excel using "table`i'.xlsx", clear allstring
	rename B run`i'
	gen nid = _n
	merge 1:1 nid using newtable, nogen
	save newtable, replace
}

drop nid

export excel "raw_one_asset_7_26_21.xlsx", replace