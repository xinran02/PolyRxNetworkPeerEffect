This repository includes R simulation code for a manuscript prepared for submission. The scripts are simulations for network peer-effect models using directional influence matrices and an undirected influence matrix.

The R scripts simulate the following scenarios:
1) When the data is generated with directional peer-effects, and the peer-effects are estimated using directional influence matrices (when the model is correctly specified);
2) When the data is generated with directional peer-effects, and a single peer-effect is estimated using an undirected influence matrix constructed from an undirected adjacency matrix (when the model is mis-specified);
3) When peer-effects only transmit through either inbound, outbound, or mutual dyads and it is how the data is generated, a single peer-effect using an undirected influence matrix is estimated.

The bias, mean-squared-error (MSE), and coverage of 95% confidence intervals are computed to assess model estimation. 
