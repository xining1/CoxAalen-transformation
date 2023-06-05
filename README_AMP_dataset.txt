

A Semiparametric Cox-Aalen Transformation Model with Censored Data

by Xi Ning, Yinghao Pan, Yanqing Sun, and Peter B.Gilbert


%%%%%%%%% Description of 'AMP_survival_20211102.csv' dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STRUCTURE
1 record per once-randomized MITT participant, i.e. our efficacy cohort as described in Corey, Gilbert et al., NEJM, 2021. 
N=4611 in HVTN 704/HPTN 085 and HVTN 703/HPTN 081 combined.

VARIABLES
-	protocol: protocol (HVTN 704 or HVTN 703)
-	pub_id: publication ID
-	tx: treatment assignment (T1, T2, or C3), 
-	tx_pool: vrc01-pooled treatment assignment (T1+T2 or C3)
-	rx: treatment label
-	hiv1event: indicator of HIV infected primary case (1=yes, 0=no)*
-	hiv1survday: primary follow-up time (days) since enrollment*
-	age: age at enrollment (years)
-	agegrp: age group at enrollment (years)
-	country: country of origin

*Follow-up time is censored at tau, the minimum time point with at least 15 participants at risk of HIV infection per treatment group and protocol, during the primary 80-week follow-up period. For HIV-infected primary cases, ‘hiv1surday’ uses the adjudicated diagnosis date of infection, which is based on ELISA. For uninfected controls, follow-up time is censored at their last negative HIV sample collection date.



%%%%%%%%%%% Instructions for using the R codes for AMP dataset %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Download the AMP dataset and open it in R. Then, install the required R packages and run the script in R.










