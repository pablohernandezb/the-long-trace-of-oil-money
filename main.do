
/* LOADING DATASETS */

import excel "C:\Users\pabherna\OneDrive - Texas Tech University\Spring 2021\POLS 5367 - International Political Economy\Research Paper\data\stata_nodes.xlsx", sheet("Sheet1") firstrow

use actors.dta, replace

/* VARIABLES */

* Generating a dummy variable for Venezuela
generate venezuela = 1 if country=="Venezuela"
replace venezuela = 0 if venezuela==.

* Generating a unified variable for cutpoint strength
gen cutpoint = 0
replace cutpoint = 1 if cutpoint_weak == 1
replace cutpoint = 2 if cutpoint_strong == 1

gen ln_total_brokerage = ln(total_brokerage)
gen ln_information = ln(information)


* Generating a factor for Centrality
factor closeness betweeness eigen totaldegree, pcf
rotate
predict centrality

* Generating a variable for the level of Brokerage
sum total_brokerage
codebook total_brokerage
hist total_brokerage
gen brokerage_level = 0
replace brokerage_level = 1 if total_brokerage>1.05
replace brokerage_level = 2 if total_brokerage>40.39

* Droping observations with no information on country
drop if transnational<0

/* VARIABLES ANALYSIS */

* Correlation between variables
pwcorr intl_trial us transnational venezuela family public oil lavajato sanctioned ///
       hr_violations cutpoint total_brokerage centrality, star(0.1)

pwcorr intl_trial us transnational venezuela family public oil lavajato sanctioned ///
       hr_violations cutpoint total_brokerage centrality, star(0.05)

pwcorr intl_trial us transnational venezuela family public oil lavajato sanctioned ///
       hr_violations cutpoint total_brokerage centrality, star(0.01)

* Summary Table for Variables
sum intl_trial us transnational venezuela family public oil lavajato sanctioned ///
    hr_violations cutpoint brokerage_level centrality

* Partial and Semipartial Correlation

pcorr intl_trial transnational venezuela family public oil lavajato sanctioned ///
      hr_violations cutpoint total_brokerage centrality

pcorr us transnational venezuela family public oil lavajato sanctioned ///
      hr_violations cutpoint total_brokerage centrality
	  
* Creating group of variables

global nationality i.transnational i.venezuela
global networks i.family i.public i.oil i.lavajato
global sanctions i.sanctioned i.hr_violations
global positioning i.cutpoint i.brokerage_level
global centrality closeness betweeness eigen indegree outdegree


/* MODELS */
	
* Model 1: Internatial Trials

* Model 1.1: Unrestricted    
logit intl_trial $nationality $networks $sanctions $positioning centrality 
estimates store model1_1

*logit intl_trial i.transnational i.venezuela i.family i.public i.oil i.sanctioned##i.brokerage_level i.lavajato i.hr_violations i.cutpoint centrality

margins, at(brokerage_level=(0 1 2))
marginsplot

margins, at(brokerage_level=(0 1 2)) by(venezuela)
marginsplot

margins, dydx(brokerage_level) by(public)
marginsplot

* Model 1.2: Nationality
logit intl_trial $nationality centrality 
estimates store model1_2

* Model 1.2: Network Ties
logit intl_trial $networks centrality 
estimates store model1_3

* Model 1.3: Sanctions 
logit intl_trial $sanctions centrality 
estimates store model1_4

* Model 1.4: Positioning 
logit intl_trial $positioning centrality 
estimates store model1_5

* Output for Logit Regressions: Long Format
esttab model1_2 model1_3 model1_4 model1_5 model1_1  using regresion_table_m1.tex, b(%9.2f) se(%9.2f) star(* 0.10 * 0.05 ** 0.01) sty(fixed) scalars(r2_p chi2 ll) label legend varlabel(_cons Constant) unstack posthead("") prefoot("") postfoot("") noomitted nobaselevels replace

* Model 2: Unrestricted Model (Trial at the US)
logit us $nationality $networks $sanctions $positioning centrality 
estimates store model2_1

margins, at(cutpoint=(0 1 2)) by(venezuela)
marginsplot

margins, at(brokerage_level=(0 1 2)) by(venezuela)
marginsplot

* Model 2.2: Nationality
logit us $nationality centrality 
estimates store model2_2

* Model 2.2: Network Ties
logit us $networks centrality 
estimates store model2_3

* Model 2.3: Sanctions 
logit us $sanctions centrality 
estimates store model2_4

* Model 2.4: Positioning 
logit us $positioning centrality 
estimates store model2_5

* Output for Logit Regressions: Long Format
esttab model2_2 model2_3 model2_4 model2_5 model2_1  using regresion_table_m2.tex, b(%9.2f) se(%9.2f) star(* 0.10 * 0.05 ** 0.01) sty(fixed) scalars(r2_p chi2 ll) label legend varlabel(_cons Constant) unstack posthead("") prefoot("") postfoot("") noomitted nobaselevels replace

* Output for Logit Regressions for Main Models: Wide Format
esttab model_1 model2_1 using regresion_table_m.tex, wide star(* 0.10 * 0.05 ** 0.01) sty(fixed) scalars(r2_p chi2 ll) label legend varlabel(_cons Constant) unstack posthead("") prefoot("") postfoot("") noomitted nobaselevels replace

*** Models on Centrality and Prestige

pwcorr centrality transnational venezuela family public oil lavajato sanctioned hr_violations

nbreg closeness $nationality $networks $sanctions
estimates store model3_1

nbreg betweeness $nationality $networks $sanctions
estimates store model3_2

nbreg totaldegree $nationality $networks $sanctions
estimates store model3_3

nbreg eigen $nationality $networks $sanctions
estimates store model3_4

reg centrality $nationality $networks $sanctions
estimates store model3_5

* Output for Centrality models: Plain Format
esttab model3_1 model3_2 model3_3 model3_5  using regresion_table_m3.tex, b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) sty(fixed) scalars() label legend varlabel(_cons Constant) unstack posthead("") prefoot("") postfoot("") noomitted nobaselevels drop(lnalpha) replace

reg closeness $nationality $networks $sanctions
estimates store model3_1

reg betweeness $nationality $networks $sanctions
estimates store model3_2

reg totaldegree $nationality $networks $sanctions
estimates store model3_3

reg eigen $nationality $networks $sanctions
estimates store model3_4

reg centrality $nationality $networks $sanctions
estimates store model3_5

* Output for Centrality models: Plain Format
esttab model3_1 model3_2 model3_3 model3_4 model3_5  using regresion_table_m3_reg.tex, b(%9.2f) se(%9.2f) star(* 0.10 ** 0.05 *** 0.01) sty(fixed) scalars() label legend varlabel(_cons Constant) unstack posthead("") prefoot("") postfoot("") noomitted nobaselevels replace


oprobit cutpoint $nationality $networks $sanctions centrality
oprobit brokerage_level $nationality $networks $sanctions centrality

*** Path Analysis

pathreg (brokerage_level transnational venezuela family public oil lavajato sanctioned hr_violations cutpoint centrality) /// 
        (cutpoint transnational venezuela family public oil lavajato sanctioned hr_violations brokerage_level centrality)  /// 
		(intl_trial transnational venezuela family public oil lavajato sanctioned hr_violations cutpoint brokerage_level centrality)

** Analysis

** Barplots with confidence intervals
* https://stats.idre.ucla.edu/stata/faq/how-can-i-make-a-bar-graph-with-error-bars/

use actors.dta, replace

keep if cutpoint!=0
tabulate intl_trial cutpoint, cell
tabulate venezuela cutpoint, cell

* Intl'Trials by Nationality and Cutpoint

use actors.dta, replace

collapse (mean) meanintl_trial=intl_trial (sd) sdintl_trial=intl_trial (count) n=intl_trial, by(cutpoint venezuela)

generate hiintl_trial = meanintl_trial + invttail(n-1,0.05)*(sdintl_trial / sqrt(n))
generate lowintl_trial = meanintl_trial - invttail(n-1,0.05)*(sdintl_trial / sqrt(n))

generate cutpointvenezuela = venezuela  if cutpoint == 0
replace cutpointvenezuela = venezuela+3  if cutpoint == 1
replace cutpointvenezuela = venezuela+6  if cutpoint == 2

sort cutpointvenezuela
list cutpointvenezuela cutpoint venezuela, sepby(cutpoint)

twoway (bar meanintl_trial cutpointvenezuela if venezuela==0) ///
       (bar meanintl_trial cutpointvenezuela if venezuela==1) ///
       (rcap hiintl_trial lowintl_trial cutpointvenezuela), ///
       legend(row(1) order(1 "International" 2 "Venezuelan") ) ///
       xlabel( 0.5 "No" 3.5 "Weak" 6.5 "Strong", noticks) ///
       xtitle("Cutpoint Level") ytitle("Mean International Trial")
	   
* With trials in the US
	   
use actors.dta, replace

collapse (mean) meanus=us (sd) sdus=us (count) n=us, by(cutpoint venezuela)

generate hius = meanus + invttail(n-1,0.1)*(sdus / sqrt(n))
generate lowus = meanus - invttail(n-1,0.1)*(sdus / sqrt(n))

generate cutpointvenezuela = venezuela  if cutpoint == 0
replace cutpointvenezuela = venezuela+3  if cutpoint == 1
replace cutpointvenezuela = venezuela+6  if cutpoint == 2

sort cutpointvenezuela
list cutpointvenezuela cutpoint venezuela, sepby(cutpoint)

twoway (bar meanus cutpointvenezuela if venezuela==0) ///
       (bar meanus cutpointvenezuela if venezuela==1) ///
       (rcap hius lowus cutpointvenezuela), ///
       legend(row(1) order(1 "International" 2 "Venezuelan") ) ///
       xlabel( 0.5 "No" 3.5 "Weak" 6.5 "Strong", noticks) ///
       xtitle("Cutpoint Level") ytitle("Mean International Trial in the US")

* Intl'Trials by Nationality and Brokerage

use actors.dta, replace

collapse (mean) meanintl_trial=intl_trial (sd) sdintl_trial=intl_trial (count) n=intl_trial, by(brokerage_level venezuela)

generate hiintl_trial = meanintl_trial + invttail(n-1,0.1)*(sdintl_trial / sqrt(n))
generate lowintl_trial = meanintl_trial - invttail(n-1,0.1)*(sdintl_trial / sqrt(n))

generate brokeragevenezuela = venezuela  if brokerage_level == 0
replace brokeragevenezuela = venezuela+3  if brokerage_level == 1
replace brokeragevenezuela = venezuela+6  if brokerage_level == 2

sort brokeragevenezuela
list brokeragevenezuela brokerage_level venezuela, sepby(brokerage_level)

twoway (bar meanintl_trial brokeragevenezuela if venezuela==0) ///
       (bar meanintl_trial brokeragevenezuela if venezuela==1) ///
       (rcap hiintl_trial lowintl_trial brokeragevenezuela), ///
       legend(row(1) order(1 "International" 2 "Venezuelan") ) ///
       xlabel( 0.5 "Low" 3.5 "Medium" 6.5 "High", noticks) ///
       xtitle("Brokerage Level") ytitle("Mean International Trial")

sort centrality	   
gen centrality_rank = _n

keep if id==1257 ///
| id==220 ///
| id==1458 ///
| id==55 ///
| id==694 ///
| id==1187 ///
| id==1318 ///
| id==50 ///
| id==1379 ///
| id==1189 ///
| id==271 ///
| id==696 ///
| id==1197 ///
| id==979 ///
| id==388 ///
| id==933 ///
| id==1543 ///
| id==80 ///
| id==1522 ///
| id==228 ///
| id==1207 ///
| id==1296 ///
| id==1165 ///
| id==140 ///
| id==24 ///
| id==1209 ///
| id==384 ///
| id==84 ///
| id==236 ///
| id==1346 ///
| id==684 ///
| id==953 ///
| id==648 ///
| id==861 ///
| id==1242 ///
| id==214 ///
| id==522 ///
| id==368 ///
| id==1180 ///
| id==955


	   
*** Saving the dataset

save "C:\Users\pabherna\OneDrive - Texas Tech University\Spring 2021\POLS 5367 - International Political Economy\Research Paper\data\actors.dta", replace


*** TO-DO LIST

* 1. Identify how is being targeted by the US
