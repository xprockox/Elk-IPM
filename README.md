# Elk-IPM
Integrated population model for elk in the northern range of Yellowstone. Model is a pre-birth-pulse three-stage matrix population model with the following stages: yearlings (0-1 year), young adults (2-13 years), and old adults (14+ years).

## Scripts

1.1.mgmt_survivalMatrices.R: Constructs survival matrices from collared elk (i=378 individuals over t=24 years). The observation matrix (y) is used in 3.5.IPM.R and consists of 1s and 0s indicating when (t) a given elk (i) was detected. A detection is either a VHF signal heard within the year, or a GPS point taken. 

1.2.mgmt_abundances.R: Constructs abundance data 
