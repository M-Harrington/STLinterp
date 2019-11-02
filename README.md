# STLinterp
A function that uses the [STLplus r package](https://cran.r-project.org/web/packages/stlplus/stlplus.pdf) that performs 
[Seasonal Trend decomposition using Loess](https://otexts.com/fpp2/stl.html) to estimate missing values.  The function 
automatically performs one of two forms of cross-validation to optimally select the tuning parameters and then returns
estimated values for the missing values.

Implicitly this method uses the assumption that the NAs in question are missing at random, and that the number of missing
values are small relative to the number of observations.

An example of how to implement the code using parallel computing is included.
