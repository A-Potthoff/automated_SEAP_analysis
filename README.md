# automated_SEAP_analysis
For the automatically generated output of SEAP abosrption analysis, we have created a short script that automatically handles the data and conducts statistical analysis.

The usage is simple:

move the sctip to the folder, where both the raw data and your design tables are found in. Open the script in an environment of your choice (ideally RStudio), adapt the values to those of your analysis and run the code. The script will create a results folder with all the results you need.

Aside the script you will also find sample data upon which this script was created.

This script:
* automatically detects the header of the excel
* recognizes the design tables using hashes the user adds to the excel
* can handle multiple treatments and automatically detects appropiate statistical tests
* is flexible with regards to the number of replicates
* automatically excludes empty wells
* works in less than 1 min.

This script was part of an assignment of the module "Biotechnology and Synthetic Biology" of the study programme of "Quantitative Biology" at the HHU Düsseldorf.
