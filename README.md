# Author
+ Cody Martin
+ University of Wisconsin-Madison
+ Department of Bacteriology
+ Anantharaman lab
+ Script created for Pensinger et al. 2022.

# Getting started
## Prepare input data
1) Data exported from plate reader software must be in comma-delimited format (`.csv`) and only consist of the raw data table with the column headers. In Excel, you should remove additional metadata and export *only* the raw data table as a `.csv` file. See `example_data.csv` provided here as the example data set used in the manuscript.
   1) Notice that the column labels are `Time` and the well coordinates. If your plate reader does not output the data in this format, it will need to be converted to this format.
   2) Additionally, notice the times are in `mm:ss` format or `hh:mm:ss` format. This is required for input but will be converted to `h` relative times, with the first time being 0.
2) A comma-delimited file for the plate setup is required. See `plate_setup.csv` provided here. The columns must be labeled first with `Row` and then the numbers 1-12, corresponding to the setup of a conventional 96-wellplate. The column beneath `Row` should be `A-H`, again corresponding to the rows of a 96-wellplate. The values in the cells are the names of the samples.
   1) Wells that should not be analyzed (either no samples were present or samples that are not part of this analysis) should be left **BLANK**.
   2) For this analysis in particular, the format of the sample name is as follows: `STRAIN_MEDIA_CONCENTRATION`. It is imperative that this format be followed here for this script to work.
   3) Be sure to name all replicates the **EXACT** same way with no additional spaces or other characters. Do not label them with incrementing numbers either.
3) Finally, ensure that the `growth_curve_statistics.py`, the data, and the plate setup files are in the same directory. This will make things easier...
## Prepare virtual environment to run script
1) Clone this repository:
   `git clone git@github.com:HryckowianLab/Pensinger_2022.git`
2) If you don't have conda, please install it from here: https://docs.conda.io/en/latest/miniconda.html. Then, create a conda environment: 
   `conda create -n platereader -c conda-forge python=3.10.2 numpy=1.22.2 pandas=1.4.1`
3) Next, activate the conda environment that was just created called `platereader`:
   `conda activate platereader`

# Running script
1) If you need help running the script, you can use the following command:
   1) `./growth_curve_statistics.py -h`
   2) This will print a help page with the command line arguments to use with the script. 
2) The bare minimum commands to run the script are as follows:
   1) `./growth_curve_statistics.py -i <data.csv> -k <plate_setup.csv> -o <output_basename>`
   2) <data.csv> = name of plate reader data file
   3) <plate_setup.csv> = name of plate reader setup file
   4) <output_basename> = the basename of the output files. All output files will be be named in the format `{output_basename}_NNN.csv` with `NNN` referring to the specific output file.

## Script arguments
1) `-i`: input data file in `.csv` format. See above for requirements.
2) `-k`: input plate reader setup file. See above for requirements.
3) `-o`: output basename that will be part of the all files generated.
4) `--noblank`: This *should* prevent background subtraction (ie blank media subtracted from the culture ODs) when processing the data. However, blanking was never implemented since growth curve shape statistics (max growth rate, lag time, doubling time, etc) should be unaffected by blanking (which is just a shift in values). Note: this does, however, mean that the reported OD values are not absolute values. Additionally, this means that this argument does nothing.
5) `--prism`: This forces the output processed data to be in a wide format that is more compatible with manual plotting in GraphPad Prism software. Do not use this command if you want to plot data with a programming language such as python or R as the processed data will be exported to a long format that is more suitable for these situations.
6) `--readme_name`: A `README.txt` file is generated every time the script is run that details the specifics of the values reported in each of the output files. This argument can be used to change the default name of this file from `README.txt` to another name if that file already exists.
7) `--ignore_treps`: Use this if you want to ignore individually processing technical replicates by themselves. **For the purposes of using this script**, the definition of technical replicates is the replicate samples within the *same* plate on the *same* day. This will result in only processing data at a biological replicate level, which I define as the average of technical replicates from a given day's experiment.

## What does the script do:
### Biological replicates
1) See above for my definition of biological replicate I use.
2) For biological replicate level processing, all (technical) replicates for each sample are summarized at every time point by the mean and standard deviation.
3) To calculate summary metrics, including final OD, max growth rate, lag time, and doubling time, the data are smoothed with a rolling average over 5 time intervals to minimize variability inherent with plate reader measurements. For the data in this manuscript, this corresponds to 150-min intervals. Here is a summary of the workflow for *each* strain:
   1) Final OD: take the average of the final 5 OD values (150 min)
   2) Max growth rate: the growth rate at every adjacent pair of time intervals $(t_1, t_2)$ is calculated as: 
      1) $\frac{\log(OD_2) - \log(OD_1)}{t_2 - t_1}$
      2) Then the max growth rate was defined as the as the max growth rate over a rolling average of 5 time intervals (150 min).
      3) The max growth rate is, therefore, calculated from an approximated slope and is not necessarily the instantaneous max growth rate.
   3) Lag time: Since the max growth rate corresponds to a time interval with the max rate of growth (and not an instantaneous time), the midpoint of the time interval with the max growth rate was used to calculate lag time: $t_{max} = \frac{t_1 + t_2}{2}$. Lag time is defined as half the time to reach the max growth rate: $t_{lag} = \frac{t_{max}}{2}$.
   4) Doubling time: $\frac{\log(2)}{M}$, where $M$ is the max growth rate as determined above.
   5) Note: in all mentions of $\log$, this corresponds to the natural log.
4) Processing the data in this biological replicate-centric manner is the default for this script and cannot be turned off.
 ### Technical replicates
 1) The above steps are performed, except treating each technical replicate as an independent unit. In other words, no summarizing the replicate samples from the same plate is performed. Otherwise, the same steps are performed *per sample per technical replicate*.
 2) This type of processing can be turned off using the `--ignore_treps` argument.

# Final details
- In the manuscript, the script was used in this format for all data processed:
- `./growth_curve_statistics.py -i <data.csv> -k <plate_setup.csv> -o <output_basename> --noblank --prism`