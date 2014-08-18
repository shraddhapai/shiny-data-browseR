
# Compute row-wise standard deviation for a matrix.

suppressMessages(library(inline))
suppressMessages(library(Rcpp))

computeCI <- function# compute 95% CI
(m) {
cppFuns <- .cpp_sd(); # hard work done in C++
getSD <- cxxfunction(signature(mymat="matrix", meanOut="numeric", sdOut="numeric"),
	body=cppFuns[["sdFun"]], plugin="Rcpp", cppFuns[["other"]])

sd <- numeric(nrow(m)); mu <- numeric(nrow(m));
print(system.time(getSD(mymat=m, sdOut=sd, meanOut=mu)))

cioffset <- (1.96*sd)/sqrt(ncol(m))

return(cbind(mu-cioffset, mu+cioffset))
}
.cpp_sd <- function() {

baseFn <- '
#include <math.h>
';

# computes HL-estimator.
mainStr <- '
NumericMatrix m(mymat);
NumericVector sd(sdOut);
NumericVector mu(meanOut);
int nrow=m.nrow(); int ncol=m.ncol();

for (unsigned int r=0; r < nrow; r++) {
	NumericVector v = m(r,_);
	// Remove NAN elements
	v.erase(std::remove_if(v.begin(), v.end(), isnan),
		v.end());
	
	double stdev; double mean;
	if (v.size()==0) {
		stdev = nan(""); mean=nan("");
	} else {
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	mean = sum/v.size();

	// difference from the mean
	std::vector<double> diff(v.size());
	std::transform(v.begin(), v.end(), diff.begin(),
               std::bind2nd(std::minus<double>(), mean));
	double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		stdev = std::sqrt(sq_sum / (v.size()-1));
	}
	sd[r] = stdev;
	mu[r] = mean;
}
';

return(list(sdFun=mainStr,other=baseFn))

}
