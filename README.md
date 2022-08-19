# KMOR_algorithm
Detecting outliers is an essential part of data analysis, and removing outliers from clusters can improve the accuracy of clustering significantly.
In this code, we implement the KMOR algorithm which simultaneously performs data clustering and outlier detection by introducing an additional "cluster" to the k-means algorithm to detect all outliers [1]. This code also runs based on an iterative method that optimizes the KMOR algorithm's objective function and establishes the iterative procedure's convergence [1]. 
In the KMOR algorithm, a data point that is at least γ×d_avg away from all the cluster centers is considered as an outlier, where γ is a multiplier and d_avg is the average distance calculated dy- namically during the clustering process [1].
## How to run
This code is implemented on two different datasets. The first dataset is called "iris" dataset and the second one is synthetically generated with 10 non_overlapping random cluster centers and gaussian distribution. 

The "synthetic_KMOR" code contains the KMOR algorithm for the synthetically generated dataset.

The "iris_KMOR" code contains the KMOR algorithm for the iris dataset.
## Results
According to the results of KMOR algorithm, Classification accuracy for the "iris" and "synthetically generated" dataset increased by 20% and 31% respectively compared to the K_means algorithm.

### Any Questions?
If you had any questions about using this code, Please contact [Sara Khalili](sarahkhalili89@gmail.com)

### Refrences
[1]	Gan G, Ng MK. K-means clustering with outlier removal. Pattern Recognition Letters. 2017 Apr 15;90:8-14.
