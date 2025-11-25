Real data analysis
=============================

Real data analysis of manufacturing email data (Section 6, Appendix E.1),
Real data analysis of conference interaction data (Appendix E.3),
Comparative STERGM analysis of manufacturing email data (Appendix E.2)

Files
-----

1. realdataFunctions.R
   - Defines helper functions for real data analysis
   - NOTE: source('real_data/realdataFunctions.R') is run in each of the following scripts

2. data_checking_man.R
   - Loads and preprocesses manufacturing email data
   - Produces Figures S2-S5 (saved to folder 'data_plots_man')

3. fit_man.R
   - Fits AR network and competing models to manufacturing email data
   - Produces Figures 1-2 (saved to folder 'fit_plots_man')
   - Produces Table 2 ('ics_man1.rds' and 'ics_man2.rds' in folder 'data')
   - NOTE: code sections 2 (model fitting) and 4 (forecast model fitting) can be skipped, results are cached

4. data_checking_sfhh.R
   - Loads and preprocesses conference interaction data
   - Produces Figures S8-S9 (saved to folder 'data_plots_sfhh')

5. fit_sfhh.R
   - Fits AR network and competing models to conference interaction data
   - Produces Figure S10 (saved to folder 'fit_plots_sfhh')
   - Produces Table S4 ('ics_sfhh.rds' in folder 'data')
   - NOTE: code section 1 (model fitting) can be skipped, results are cached

6. fit_man_tergm.R
   - Fits comparative STERGM models to manufacturing email data
   - Produces Figures S6-7 (saved to folder 'fit_plots_man')

Usage
-----

2-6.
  - all scripts run independently inside 'arnetworks.RProj'
  - all scripts source helper functions in 'realdataFunctions.R'
  - results from code chunks (3-2, 3-4, 6-1) which take over 5 minutes to run have been cached in folder 'data' and can be skipped
  - manufacturing email data is imported from R package networkDynamicData
  - conference interaction data is included in the repository ('tij_SFHH.dat_.gz')
  - Figures in the manuscript are saved in individual folders (see above)
  - Tables in the manuscript are saved in folder 'data'

Requirements
------------

- R 4.2.1
- Packages: arnetworks, latex2exp, networkDynamicData, pROC, statnet
