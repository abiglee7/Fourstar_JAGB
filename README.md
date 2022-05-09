# Fourstar_JAGB

Author: Abigail J. Lee

These scripts and jupyter notebooks clean, calibrate, and merge raw photometry files as created from DAOPHOT. prep_photometry.py contains all the necessary python code and the Jupyter notebooks contain examples for all 13 galaxies. 

1. clean_galaxies.ipynb cleans the photometry based off J-band chi, sharp, and error values. 
2. double check 2MASS zeropoints.ipynb then matches the brighest stars in the photometry files to 2MASS to see if there are any additional photometric offseets. Andy did prelimary matchings but I found some additional offsets on the order of <0.05 mag for several epochs. 
3. Merge Epochs.ipynb computes weighted average photometry across all individual epochs. 
