# QuantifyingWaveAttenuation
Data and scripts to accompany Pinsky et al. 2013 Ecosphere doi: 10.1890/ES13-00080.1

This repository was started in September 2022 in the interests of open science and to facilitate further research using these data. Because this is roughly a decade after the analysis was completed, the notes are starting out rather thin and there are likely to be mistakes. Files were named with dates since this was originally made without version control. Feel free to be in touch if you have questions or find errors.

Malin Pinsky
malin.pinsky@rutgers.edu

## data
- CoastProt Quant 110703.csv: The main table of data extracted from studies.
- SppData 110829.csv: data on species.

## scripts
- CoastProt Quant 130612.r: 
  - Preps data for the Matlab scripts. Reads data/CoastProt Quant 110703.csv and output/SppSummary_2011-07-04.csv. Writes output/quant_DATE.csv. 
  - Calculates sample size from output/quant_DATE.csv.
  - Finds missing data, writes output/Missing_dataDATE.csv
  - Reads in Matlab output from analysis/Fit_Cd 110622/VegData2Out_25-Mar-2012.csv. Summarizes, cleans up, and writes output/quantcd_allHDATE.csv (full data file). Also writes output/quantcdlf_allHDATE.csv after aggregating by LF-hab3-study. Writes output/quantcdhab2_allHDATE.csv aggregated by hab2-study. Writes output/quantcdhab3_allHDATE.csv aggregated to hab3-study.
  - Outputs a table of the studies by reading output/quantcd_allH2012-03-25.csv and then writing output/StudiesDATE.csv
  - Makes some of the figures for the paper. 
- Matlab scripts produce analysis/Fit_Cd 110622/VegData2Out_25-Mar-2012.csv

## output

## figures