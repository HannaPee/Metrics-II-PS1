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





































