This is file contain instructions on how to run the experiments in the paper:
"Sequential nonparametric tests for a change in distribution: an application to detecting radiological anomalies".

- "Paper_Table1.R"
This file allows to get the results in Table 1 in the paper.  The average delay time  for method X is saved to array:
"average_delay_time_X".  From these arrays one can manually get the results in Table 1.

-For  replicating Table 2 one can first run "Paper_generate_radition_ground_truth.R”, here one needs to specify the parameters of the anomaly. Then run: “Paper_simulated_radition.R”  to obtain comparisons between different methods.


-"Paper_small_anomaly.R"  is the code that can be used for obtaining the results in Table 3 in the paper. 
