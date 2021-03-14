# cs129finalProject
Final project for the CS129 class
<br />
Main fitting code is located here:<br />
./mainFittingCode/FitAllPowers.m<br />
Running it will load 38 different optical transimission and SHG data and apply the nonlinear regression model to the data to learn parameters. Output data plots during my run were automatically saved in the following folder:<br />
./mainFittingCode/FitLineshapes20210313T191602<br />
Data was exported here:<br />
./mainFittingCode/FitLineshapes20210313T191602/AllData2<br />
<br /><br />

For further data processing and error evaluation I use the following folder:<br />
./nonlinearRegressionTests<br />
I use the output of the previous section, local copy located here:<br />
./nonlinearRegressionTests/Data<br />
Script for final figure generation:<br />
./nonlinearRegressionTests/getProjectData.m<br />
Normalized RMSE calculation:<br />
./nonlinearRegressionTests/calculateNRMSEscirpt.m<br />
Results exported to:<br />
./nonlinearRegressionTests/final_plots<br />
<br /><br />
./exampleFits contains example nonlinear fits to the experimental data
