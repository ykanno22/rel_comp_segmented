Copyright (C) Yoshihiro Kanno, 2023. All rights reserved.
===============================================================================
Y. Kanno: "Confidence bound for structural response using segmented least squares: a mixed-integer programming approach." Submitted to Japan Journal of Industrial and Applied Mathematics.

The codes can run on MATLAB ver.23.2 with CPLEX ver.12.9. 

To reproduce the result presented in section 6.2, run the following codes sequentially: 
  - tri_modulus_1d_data_set.m
  - segmented_regress.m
Then run the following codes:
  - displ_analysis_loop.m
