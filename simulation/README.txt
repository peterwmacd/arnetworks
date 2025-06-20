Transitivity Model Simulation
=============================

Estimation and inference for the transitivity model.

Files
-----

1. Step1_estTransitivity.R  (main step)
   - Runs simulations from index start_try to end_try.  
   - Saves estimation and inference results as out_<tag>_<start>_<end>.rds

2. run_simjob.sbatch  
   - SLURM script to run Step1_estTransitivity.R in parallel.  
   - Usage: sbatch run_simjob.sbatch <tag>

3. Step2_combine_result.R  
   - Combines all partial results into one file: out_all_<tag>.rds

4. Step3_summary.R  
   - Summarizes the results (e.g., errors, coverage)

Usage
-----

1. Load input RData file under /Setting with correct tag.  
2. Run: sbatch run_simjob.sbatch <tag>  
3. After jobs finish, run Step2_combine_result.R  
4. Then run Step3_summary.R for analysis

Requirements
------------

- R 4.2.1  
- Packages: arnetworks, doParallel, foreach
