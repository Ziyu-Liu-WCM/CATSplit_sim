# CATSplit Simulation

The latest code is located at the **poisson** folder. All you need to get the simulation working are the Load_Simulator.R and Simulation_Main.R. The Load_Simulator.R file is to generate simulation data and save it at the **input** folder. The Simulation_Main.R is used for evaluating performance. For simple testing, you can ignore the code for command line input and simply de-comment the **For Testing** part and modify them with your desired input.

I think for now you should be super vigilant about the nPerm parameter because I haven't incorporated it in the R^2 storage process, so every time when you want to test a different nPerm for a same dataset you need to delete all results stored.

dataSet can be one of "SCI", "IBD"

featureSet can be one of "lowCor", "highCor", "leastAbun", "mostAbun", "leastVar", "mostVar"

metric can be one of "euclidean", "Weighted UniFrac", "Unweighted UniFrac", "robust"

