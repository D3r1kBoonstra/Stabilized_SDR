1. Implementation.

(1a) "~/EnvelopeAlgorithms" is collection of MATLAB functions for estimating an envelope for various multivariate statistical problems. This folder contains the computation algorithms for obtaining a generic envelope. Specifically, the subfolder "/EnvelopeAlgorithms/4_Gmanifold" contains the functions (objective functions, derivatives, etc.) for optimizing the ENDS-type optimization.

(1b) "/ENDS" contains various MATLAB functions for the ENDS method such as prediction, classification, dimension determination, subspace estimation, cross-validation tuning parameter selection, etc.



2. Simulation Examples in the article.
We use the simulation example (Q2) in the paper to demonstrate classification and prediction, dimension selection, and subspace estimation using the proposed ENDS method. Results in Table 1, 2 and 3 can be reproduced by these codes.

(2a) "Example_Classification_Prediction.m"
Includes comparisons of ENDS estimators (ENDS-LDA, ENDS-QDA) with SIR-LDA, SIR-QDA, SAVE-LDA, SAVE-QDA, DR-LDA, DR-QDA, Naive Bayes, and Support Vector Machine. We also included the Bayes error.

(2b) "Example_Dimension_Selection.m"
Demonstrates the dimension selection of ENDS using AIC and BIC.

(2c) "Example_Subspace_Estimation.m"
Includes comparisons of ENDS versus SIR, SAVE and DR in terms of subspace estimation error.