# Elk-IPM
Integrated population model for elk in the northern range of Yellowstone. Model is a pre-birth-pulse, female-only, three-stage matrix population model with the following stages: yearlings (0-1 year), young adults (2-13 years), and old adults (14+ years). A "year" as defined in this model runs from June 1 - May 31 of the following year. 

<img loading="lazy" width="1000px" src="./figures/Screenshot 2025-11-14 at 11.47.50.png" alt="image_name png" />

## The model

The full elk IPM is contained in 3.5.IPM.R and comprises four sub-models: (1) fecundity, (2) calf survival, (3) adult survival, and (4) growth from young adult to adult. These all estimate parameters to be used in the full state-space model which estimates stage-specific abundances.

<img loading="lazy" width="1000px" src="./figures/Screenshot 2025-11-14 at 11.44.55.png" alt="image_name png" />

## Data

**A full description of all data contained in this repository can be found in data/data_sources_key.txt**

The elk IPM incorporates data from hunter harvests (occuring Oct - Dec of a given year), collared elk (collared in Dec - Jan of a given year), total elk counts (estimated from aerial flights conducted in Jan - Feb), and elk classification surveys (also from aerial counts, but conducted in Mar - Apr). Hunter harvest data includes the age structure and pregnancy status of all harvested cows reported for GMUs 313 and 316 (immediately north of the Yellowstone park boundary). Collared individuals are used in the survival model (Cormack Jolly Seber), and the pregnancy status + age of each collared elk is also used. Classification surveys provide information on the percentage of cows seen in various groups (population-level average % cows) and calf:cow ratios (also population-level averages). 

<img loading="lazy" width="1000px" src="./figures/Screenshot 2025-11-14 at 11.49.05.png" alt="image_name png" />

## Data management scripts

1.1.mgmt_survivalMatrices.R: Constructs survival matrices from collared elk (i=378 individuals over t=24 years). The observation matrix (y) is used in 3.5.IPM.R and consists of 1s and 0s indicating when (t) a given elk (i) was detected. A detection is either a VHF signal heard within the year, or a GPS point taken. 
+  adultSurival_cjsMatrices.rData: the entire R environment created by "1.1.mgmt_survivalMatrices.R" including:
	- is_class1: a matrix of 378 collared individuals over 26 years. 1s represent when a given individual was in the first young adult class (aged 2-13) and 0s represent when a given individual was NOT in the first young adult class.
	- is__class2: a matrix of 378 collared individuals over 26 years. 1s represent when a given individual was in the second young adult class (aged 14+) and 0s represent when a given individual was NOT in the second young adult class.
	- y: the observation matrix of our 378 collared individuals over 26 years. 1s represent when an individual was "observed" (i.e. a GPS point was taken or a VHF signal was heard at some point during that year). 0s represent when an individual was not observed, either because they had not yet been collared, had a collar failure, or died.
	- z: the latent (true) state matrix of our 378 collared individuals over 26 years. This is essentially just y, but with additional 1s where the individual was known to be alive (e.g., in the case that a collar failed for only one year but then came back online, we would know the individual was alive during that year the collar failed). This matrix actually isn't used in the analyses, but it could be useful in the future.
	- first_seen: a vector of 378 length indicating the year (numbered 1-26) that each individual was first seen (used for censoring individuals in the CJS model prior to their capture).

1.2.mgmt_abundances.R: Constructs abundance data as follows:
+  abundanceEstimates_stages.csv: Created by the script "1.2.mgmt_abundances.R". Columns are:
	- X: unique ID for each year (artifact of the way the data were saved from the R script, not important)
	- year: year
	- n_total, n_total_LCI, n_total_UCI: the abundance estimate from flight surveys with lower and upper 95% CIs
	- percent_cow, percent_calf, ... percent_bull: percent of each demographic class from classification flights 
	- n_cow, n_calf, ... n_bull: an estimate of the number of cows, calves, etc. calculated as n_total * the percentage from the classification data
	- sigma_tot_log: standard deviation of the flight estimate (n_total) derived from n_total_L95 and n_total_U95 as: 
	sigma_tot_log = round((log(n_total_U95) - log(n_total_L95)) / (2*1.96),2))
	- prop_old_ad, prop_young_ad: from the hunter harvest age distribution data, the proportion of young adults that are either old (14+ y.o.) or young (2-13 y.o.)
	- n_cow_youngadult, n_cow_oldadult: an estimate of the number of young vs. old cows using previous "n_cow" and "prop_old_ad" or "prop_young_ad"

1.2.mgmt_fecundity.R: Constructs fecundity data as follows:
+  fecundity.csv: Created by the script "1.3.mgmt_fecundity.R". Columns are:
	- X: unique ID for each year (artifact of the way the data were saved from the R script, not important)
	- year: year
	- n_total, n_total_LCI, n_total_UCI: the abundance estimate from flight surveys with lower and upper 95% CIs
	- calf_cow_ratio, spike_cow_ratio ... bull_cow_ratio: the number of calves, bulls, etc. for every cow (based on classification flights)
	- percent_cows: the percent of the population that are cows (based on classification flights)
	- n_cow, n_calf, ... n_bull: an estimate of the number of cows, calves, etc. calculated as n_total * the percentage from the classification data
	- percent_cows_old, percent_cows_yearling, percent_cows_young: using the hunter harvest age distribution data, these are estimates of the proportion of harvest elk in each age class
	- yearling_num_preg, young_num_preg, old_num_preg: the number of each age class that tested positive (PSPB+) during capture. For years where captures did not take place, this is estimated from hunter reports from GMUs directly north of the park
	- yearling_num_capt, young_num_capt, old_num_capt: the number of each age class that were captured each year. For years where captures did not take place, this number represents the number of cows in each stage that were harvested in GMUs directly north of the park
	- yearling_prop_pregnant, young_prop_pregnant, old_prop_pregnant: number pregnant / number captured from previous columns
	- n_cows_young, n_cows_old: an estimate of the number of young vs. old cows using previous "n_cow" and "percent_cows_old", "percent_cows_young"
	- expected_calves_from_young, expected_calves_from_old: proportion pregnant for each stage * n_cows in each stage (assumes no abortion--if a cow is pregnant, they give birth)
	- expected_calves: sum of expected_calves_from_young + expected_calves_from_old
	- calf_surv: the number of calves counted during classification surveys (n_calf) divided by expected_calves
	- harvested_age13: the proportion of hunter-harvested individuals that were 13 years old
	- harvested_total: the total number of young adults harvested
	- prop_age13: harvested_age13 divided by harvested_total
	- n_age13: prop_age13 * n_cows_young
