# StabilityMeasures


Source code for manuscript "Selection of variables for multivariable models: opportunities and limitations in quantifying model stability by resampling"
by Christine Wallisch, Daniela Dunkler, Geraldine Rauch, Riccardo de Bin, Georg Heinze, Statistics in Medicine 2020, DOI:10.1002/sim.8779.

For questions, comments or remarks about the code please contact Christine Wallisch (christine.wallisch@meduniwien.ac.at).


These are the proceeding steps to reproduce the results presented in the article:

Note: The simulation study was conducted on the Vienna Scientific Cluster. The code attached here has been modified to run on a single PC. Thus, results may be slightly different.

## Simulation study

### LINEAR REGRESSION 

See folder 'simulation study/1 code for linear regression'

1. Run 'sim_setup.R' to save a matrix with all simulation szenarios
2. Run 'sim_results_beta.R' to obtain all coefficients in the simulation szenarios. (This takes several weeks on a single PC.)
3. Run 'sim_measures.R' to calculate stability measures from the coefficients. 
4. Run 'Table2.R', 'Table3.R', 'Figure1.R', 'Figure2.R' to obtain the results presented in the article.
5. Run 'Table2.R', 'Supplementary_material_S3.Rmd', 'Supplementary_material_S4-5.R', 'Table 3.R', 'Supplementary_material_S6-7.R' to obtain the results in Supplementary Material S2-S7.

'function_setup.R', 'function_beta.R', 'function_glmnetAIC.R' and 'function_measures.R' are auxiliary files.


### LOGISTIC REGRESSION

Repeat the first 3 steps described above but use the files in folder '2 code for logistic regression'.
Then, run 'Supplementary_Figure7.R', 'Supplementary_Table3.R', 'Supplementary_Figure8.R' to obtain the results in Supplementary Material S8.


### COX REGRESSION

Repeat the first 3 steps described above but use the files in folder '3 code for Cox regression'.
Then, run 'Supplementary_Figure9.R', 'Supplementary_Table4.R', 'Supplementary_Figure10.R' to obtain the results in Supplementary Material S9.

## Example

The code for the EXAMPLE on CARDIOVASCULAR DISEASE is given in the folder 'example CVD'. 
As the original data set is sensitive and cannot be published, we synthesized the data set with the R package 'synthpop'. 
Therefore, the results with the synthesized data set will not match the results presented in the article that were obtained with the original data set.
Run 'modelling_and_measures.R' to obtain reproduce the analysis of Table 4 in the article (which is based on the original data), and to reproduce Supplementary Table 5 for the synthesized data set.
