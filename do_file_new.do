/// Metrics II PS 1 /// 

** housekeeping
clear all                   // remove anything old stored
set more off, permanently   // tell Stata not to pause
set linesize 255            // set line length for the log file
version                     // check the version of the command interpreter

* Set working directory to the current repo folder
cd "C:\Users\42610\OneDrive - Handelshögskolan i Stockholm\Documents\Metrics_II_PS1"
global wd "`c(pwd)'"

* Create folders if they do not exist
cap mkdir figures
cap mkdir output
cap mkdir logs

** Create a RED.ME file, and choose the name (in this case test_stata)

! echo #  metrics_ii_ps1 >> README.md

**Initialize git code 

! git init

** Add READ.ME file and comit 

! git add README.md
! git commit -m "first upload"

! git branch -M main
! git push -u origin main

** Define the directory where we wanto add this file 

! git remote add origin https://github.com/HannaPee/Metrics-II-PS1.git

** capture
cap log close // close a log-file, if one is open
log using "metrics_ii_ps1.log", replace


****** Question 1 **********

** simulate data 

set obs 100 
set seed 2604
gen i = _n
gen X = 1 if i <= 50
	replace X = 2 if i > 50 & i <= 75 
	replace X = 3 if 75 < i

	
* simulate the treatment variable * 	
gen u = runiform()
gen D = u <= 1/3 if X == 1 
replace D = u <= 2/3 if X == 2 
replace D = 1 if X == 3 
drop u 

gen y_0= runiform()
gen y_1= runiform() * X


*** generate indiivdual causal effects ***
gen tau_i=y_1 - y_0

*** plot the individual causak effects *** 

twoway ///
    (kdensity tau_i if X==1, lcolor(navy)   lwidth(medthick)) ///
    (kdensity tau_i if X==2, lcolor(maroon) lwidth(medthick)) ///
    (kdensity tau_i if X==3, lcolor(forest_green) lwidth(medthick)), ///
    legend(order(1 "X = 1" 2 "X = 2" 3 "X = 3")) ///
    xtitle("tau") ///
	ytitle("Kernel density")
 

graph export "figures/kdensity_tau.pdf", replace


** Estimation ** 

gen y = cond(D == 1, y_1, y_0)

eststo clear
eststo att: reg y D
forvalues i = 1/3 {
    eststo catt_`i': reg y D if X == `i'
}

esttab att catt_1 catt_2 catt_3 using estimation.tex, ///
    nostar mtitles("ATE/ATT" "X=1" "X=2" "X=3") replace

** Calculate the true effects ** 
qui sum tau_i
scalar ate = r(mean)

qui sum tau_i if D==1
scalar att = r(mean)

forvalues i = 1/3 {
    qui sum tau_i if D==1 & X==`i'
    scalar catt_`i' = r(mean)
}

eststo clear
eststo att_table
estadd scalar ATE = ate
estadd scalar ATT = att
estadd scalar CATT1 = catt_1
estadd scalar CATT2 = catt_2
estadd scalar CATT3 = catt_3

** Export table ** 
esttab att_table using att_table.tex, ///
    cells(none) ///
    stats(ATE ATT CATT1 CATT2 CATT3, ///
          labels("ATE" "ATT" "CATT1" "CATT2" "CATT3")) ///
    replace



! git add figures/
! git commit -m "Add figures"
! git push


********** Question 2 *********


*generate data * 

set obs 200 
set seed 2604
gen Y_0 = rnormal()
gen X = rbinomial(1,0.3)

local a0 = 1
local a1 = -2

gen tau_i = `a0' + `a1' * X

gen u = uniform()
gen D = u > 0.5 
drop u 

gen Y = Y_0 + tau_i*D

* Estimate tau_hat by comparing means*

estpost summarize Y if D==0
eststo m0

estpost summarize Y if D==1
eststo m1

esttab m0 m1 using mean_table.tex, ///
	cells("mean") ///
    mtitles("D=0" "D=1") ///
    replace ///
    booktabs

* Difference in means via regression
reg Y D
eststo diff

* Export table
esttab diff using mean_2_table.tex, ///
	nostar mtitles("ATT/ATE") replace

*** Estimate two CATT ** 
preserve

collapse (mean) Y, by(X D)

reshape wide Y, i(X) j(D)

gen catt = Y1 - Y0

list catt if X == 0 
list catt if X == 1 
	
restore	

* Difference in means via regression
reg Y D if X==0
eststo catt0

reg Y D if X==1
eststo catt1

esttab catt0 catt1 using catt_reg_table.tex, ///
    cells("b") ///
    keep(D) ///
    mtitles("X=0" "X=1") ///
    replace ///
    booktabs

	
****** Questoin 3 ******

clear all

set obs 100 
set seed 2604

gen Y_1 = rnormal(2,1)
gen Y_0= rnormal(1,1)
gen tau_i = Y_1 - Y_0
gen D = runiform() > 0.5 
gen Y= Y_1*D + Y_0*(1-D)

* plot the taus* 

histogram tau_i, ///
    bin(20) ///
    xtitle("tau_i") ///
    ytitle("Frequency")

 graph export "figures/tau_histogram.pdf", replace
 
 *drop all variables except Y and D* 
 drop Y_1
 drop Y_0 
 drop tau_i 
 
 *regress and save Y on D* 
 
 reg Y D 
 scalar tau_true = _b[D] 

 **** randomized inference****
 
 * generate potential outcomes once*
 
gen Y_0_h0= rnormal(1,1)
gen Y_1_h0 = Y_0_h0

* Create a frame to store the simulated coefficients
frame create results tau

* generating variables
gen D_a = .
gen Y_a = .

*looping over randomdized assignments*

forvalues j = 1/1000 {
    replace D_a = runiform() > 0.5
    replace Y_a = Y_0_h0 + (Y_1_h0 - Y_0_h0)*D_a

    quietly reg Y_a D_a
    local tauhat = _b[D_a]
    frame post results (`tauhat')
}

* Work with the results frame
frame change results

histogram tau, ///
    bin(30) ///
    xline(`=tau_true') ///
    xtitle("Tau hat") ///
    ytitle("Frequency") ///

graph export "figures/randomization_distribution.pdf", replace

* Count if taus are larger than true tau* 

count if tau > tau_true
display r(N)/_N


**** bootstrap ****
frame change default

frame create bootsrp tau_boot

qui forval i = 1/1000 {
	preserve 
		bsample 
		
		qui reg Y D

		local tauhat = _b[D]
		frame post bootsrp (`tauhat')
		
	restore
}

frame change bootsrp

* plot the bootstrap taus* 

histogram tau_boot, ///
    bin(30) ///
    xline(`=tau_true') ///
    xtitle("Tau_hat_bootstrap") ///
    ytitle("Frequency")
	
graph export "figures/bootstrap_distribution.pdf", replace

* Count if bootstrap taus are larger than true tau* 

count if tau_boot < 0
display r(N)/_N

! git add do_file [Recovered).do














