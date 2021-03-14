# cs129finalProject
Final project for the CS129 class

Main fitting code is located here:
./mainFittingCode/FitAllPowers.m
Running it will load 38 different optical transimission and SHG data and apply the nonlinear regression model to the data to learn parameters. Output data plots during my run were automatically saved in the following folder:
./mainFittingCode/FitLineshapes20210313T191602
Data was exported here:
./mainFittingCode/FitLineshapes20210313T191602/AllData2


For further data processing and error evaluation I use the following folder:
./nonlinearRegressionTests
I use the output of the previous section, local copy located here:
./nonlinearRegressionTests/Data
Script for final figure generation:
./nonlinearRegressionTests/getProjectData.m
Normalized RMSE calculation:
./nonlinearRegressionTests/calculateNRMSEscirpt.m
Results exported to:
./nonlinearRegressionTests/final_plots

./exampleFits contains example nonlinear fits to the experimental data
