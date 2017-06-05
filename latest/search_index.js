var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "(Image: GeoStatsLogo)(Image: Build Status) (Image: GeoStats) (Image: Coverage Status)"
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "This package provides high-performance implementations of geostatistical algorithms for the Julia programming language. It is in its initial development, and currently only implements Kriging estimation methods. More features will be added as the Julia type system matures."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Get the latest stable release with Julia's package manager:Pkg.add(\"GeoStats\")"
},

{
    "location": "index.html#Quick-example-1",
    "page": "Home",
    "title": "Quick example",
    "category": "section",
    "text": "Below is a quick example of usage:using GeoStats\nsrand(2017) # hide\n\n# create some data\ndim, nobs = 3, 10\nX = rand(dim, nobs); z = rand(nobs)\n\n# target location\nxₒ = rand(dim)\n\n# define a variogram model\nγ = GaussianVariogram(1.,1.,0.) # sill, range and nugget\n\n# define an estimator (i.e. build the Kriging system)\nsimkrig = SimpleKriging(X, z, γ, mean(z))\nordkrig = OrdinaryKriging(X, z, γ)\nunikrig = UniversalKriging(X, z, γ, 1)\n\n# estimate at target location\nμ, σ² = estimate(simkrig, xₒ)\nprintln(\"Simple Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide\nμ, σ² = estimate(ordkrig, xₒ)\nprintln(\"Ordinary Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide\nμ, σ² = estimate(unikrig, xₒ)\nprintln(\"Universal Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide"
},

{
    "location": "variograms.html#",
    "page": "Variograms",
    "title": "Variograms",
    "category": "page",
    "text": ""
},

{
    "location": "variograms.html#Variograms-1",
    "page": "Variograms",
    "title": "Variograms",
    "category": "section",
    "text": "newcommandxboldsymbolx\nnewcommand1mathbb1In an intrinsically stationary isotropic model, the variogram is only a function of the distance between any two points x_1x_2 in mathbbR^m:gamma(x_1x_2) = gamma(x_1 - x_2) = gamma(h)The same holds for the covariance, which is directly related gamma(h) = cov(0) - cov(h). This package implements a few commonly used stationary models:"
},

{
    "location": "variograms.html#GeoStats.GaussianVariogram",
    "page": "Variograms",
    "title": "GeoStats.GaussianVariogram",
    "category": "Type",
    "text": "GaussianVariogram(s, r, n)\n\nINPUTS:\n\ns ∈ ℜ - sill\nr ∈ ℜ - range\nn ∈ ℜ - nugget\n\n\n\n"
},

{
    "location": "variograms.html#Gaussian-1",
    "page": "Variograms",
    "title": "Gaussian",
    "category": "section",
    "text": "gamma(h) = (s - n) left1 - expleft(-left(frachrright)^2right)right + n cdot 1_(0infty)(h)GaussianVariogram"
},

{
    "location": "variograms.html#GeoStats.SphericalVariogram",
    "page": "Variograms",
    "title": "GeoStats.SphericalVariogram",
    "category": "Type",
    "text": "SphericalVariogram(s, r, n)\n\nINPUTS:\n\ns ∈ ℜ - sill\nr ∈ ℜ - range\nn ∈ ℜ - nugget\n\n\n\n"
},

{
    "location": "variograms.html#Spherical-1",
    "page": "Variograms",
    "title": "Spherical",
    "category": "section",
    "text": "gamma(h) = (s - n) leftleft(frac32left(frachrright) + frac12left(frachrright)^3right) cdot 1_(0r)(h) + 1_rinfty)(h)right + n cdot 1_(0infty)(h)SphericalVariogram"
},

{
    "location": "variograms.html#GeoStats.ExponentialVariogram",
    "page": "Variograms",
    "title": "GeoStats.ExponentialVariogram",
    "category": "Type",
    "text": "ExponentialVariogram(s, r, n)\n\nINPUTS:\n\ns ∈ ℜ - sill\nr ∈ ℜ - range\nn ∈ ℜ - nugget\n\n\n\n"
},

{
    "location": "variograms.html#Exponential-1",
    "page": "Variograms",
    "title": "Exponential",
    "category": "section",
    "text": "gamma(h) = (s - n) left1 - expleft(-frachrright)right + n cdot 1_(0infty)(h)\nExponentialVariogram"
},

{
    "location": "estimation.html#",
    "page": "Estimation",
    "title": "Estimation",
    "category": "page",
    "text": ""
},

{
    "location": "estimation.html#GeoStats.estimate",
    "page": "Estimation",
    "title": "GeoStats.estimate",
    "category": "Function",
    "text": "estimate(estimator, xₒ)\n\nCompute mean and variance for the estimator at location xₒ.\n\n\n\n"
},

{
    "location": "estimation.html#GeoStats.fit!",
    "page": "Estimation",
    "title": "GeoStats.fit!",
    "category": "Function",
    "text": "fit!(estimator, X, z)\n\nBuild Kriging system from locations X with values z and save factorization in estimator.\n\n\n\n"
},

{
    "location": "estimation.html#GeoStats.weights",
    "page": "Estimation",
    "title": "GeoStats.weights",
    "category": "Function",
    "text": "weights(estimator, xₒ)\n\nCompute the weights λ (and Lagrange multipliers ν) for the estimator at location xₒ.\n\n\n\n"
},

{
    "location": "estimation.html#Estimation-1",
    "page": "Estimation",
    "title": "Estimation",
    "category": "section",
    "text": "A Kriging estimator has the form:newcommandxboldsymbolx\nhatZ(x_0) = lambda_1 Z(x_1) + lambda_2 Z(x_2) + cdots + lambda_n Z(x_n)quad x_i in mathbbR^m lambda_i in mathbbRwith Zcolon mathbbR^m times Omega to mathbbR a random field.This package implements the following Kriging variants:Simple Kriging\nOrdinary Kriging\nUniversal Kriging (polynomial drift for the mean)All these variants follow the same interface: an estimator object is first created with a given data configuration and variogram model, and then estimates are made at various locations.The object construction takes care of building the Kriging system and factorizing the LHS with an appropriate decomposition (e.g. Cholesky, LU). The estimate method performs the estimation at a given location:estimateA typical use of the interface is as follows:# build and factorize the system\nsimkrig = SimpleKriging(X, z, cov, mean(z))\n\n# estimate at various locations\nfor xₒ in locations\n  μ, σ² = estimate(simkrig, xₒ)\nendIn case the data configuration needs to be changed in a loop (e.g. sequential Gaussian simulation), one can keep all the parameters fixed and only update the factorization with the fit! method:fit!For advanced users, the Kriging weights and Lagrange multipliers at a given location can be accessed with the weights method. This method returns an AbstractWeights object containing a field λ for the weights and a field ν for the Lagrange multipliers:weights# weights and Lagrange multipliers\nOKweights = weights(ordkrig, xₒ)\nOKweights.λ, OKweights.ν"
},

{
    "location": "estimation.html#GeoStats.SimpleKriging",
    "page": "Estimation",
    "title": "GeoStats.SimpleKriging",
    "category": "Type",
    "text": "SimpleKriging(X, z, γ, μ)\n\nINPUTS:\n\nX ∈ ℜ^(mxn) - matrix of data locations\nz ∈ ℜⁿ      - vector of observations for X\nγ           - variogram model\nμ ∈ ℜ       - mean of z\n\n\n\n"
},

{
    "location": "estimation.html#Simple-Kriging-1",
    "page": "Estimation",
    "title": "Simple Kriging",
    "category": "section",
    "text": "In Simple Kriging, the mean mu of the random field is assumed to be constant and known. The resulting linear system is:newcommandCboldsymbolC\nnewcommandcboldsymbolc\nnewcommandlboldsymbollambda\nnewcommand1boldsymbol1\nnewcommandzboldsymbolz\nbeginbmatrix\ncov(x_1x_2)  cov(x_1x_2)  cdots  cov(x_1x_n) \ncov(x_2x_1)  cov(x_2x_2)  cdots  cov(x_2x_n) \nvdots  vdots  ddots  vdots \ncov(x_nx_1)  cov(x_nx_2)  cdots  cov(x_nx_n)\nendbmatrix\nbeginbmatrix\nlambda_1 \nlambda_2 \nvdots \nlambda_n\nendbmatrix\n=\nbeginbmatrix\ncov(x_1x_0) \ncov(x_2x_0) \nvdots \ncov(x_nx_0)\nendbmatrixor in matricial form Cl = c. We subtract the given mean from the observations boldsymboly = z - mu 1 and compute the mean and variance at location x_0:mu(x_0) = mu + boldsymboly^top lsigma^2(x_0) = cov(x_0x_0) - c^top lSimpleKriging"
},

{
    "location": "estimation.html#GeoStats.OrdinaryKriging",
    "page": "Estimation",
    "title": "GeoStats.OrdinaryKriging",
    "category": "Type",
    "text": "OrdinaryKriging(X, z, γ)\n\nINPUTS:\n\nX ∈ ℜ^(mxn) - matrix of data locations\nz ∈ ℜⁿ      - vector of observations for X\nγ           - variogram model\n\n\n\n"
},

{
    "location": "estimation.html#Ordinary-Kriging-1",
    "page": "Estimation",
    "title": "Ordinary Kriging",
    "category": "section",
    "text": "In Ordinary Kriging the mean of the random field is assumed to be constant and unknown. The resulting linear system is:beginbmatrix\nC  1 \n1^top  0\nendbmatrix\nbeginbmatrix\nl \nnu\nendbmatrix\n=\nbeginbmatrix\nc \n1\nendbmatrixwith nu the Lagrange multiplier associated with the constraint 1^top l = 1. The mean and variance at location x_0 are given by:mu(x_0) = z^top lambdasigma^2(x_0) =  cov(x_0x_0) - beginbmatrix c  1 endbmatrix^top beginbmatrix l  nu endbmatrixOrdinaryKriging"
},

{
    "location": "estimation.html#GeoStats.UniversalKriging",
    "page": "Estimation",
    "title": "GeoStats.UniversalKriging",
    "category": "Type",
    "text": "UniversalKriging(X, z, γ, degree)\n\nINPUTS:\n\nX ∈ ℜ^(mxn) - matrix of data locations\nz ∈ ℜⁿ      - vector of observations for X\nγ           - variogram model\ndegree      - polynomial degree for the mean\n\nOrdinary Kriging is recovered for 0th degree polynomial.\n\n\n\n"
},

{
    "location": "estimation.html#Universal-Kriging-1",
    "page": "Estimation",
    "title": "Universal Kriging",
    "category": "section",
    "text": "In Universal Kriging, the mean of the random field is assumed to be a polynomial:mu(x) = sum_k=1^N_d beta_k f_k(x)with N_d monomials f_k of degree up to d. For example, in 2D there are 6 monomials of degree up to 2:mu(x_1x_2) =  beta_1 1 + beta_2 x_1 + beta_3 x_2 + beta_4 x_1 x_2 + beta_5 x_1^2 + beta_6 x_2^2The choice of the degree d determines the size of the polynomial matrixnewcommandFboldsymbolF\nnewcommandfboldsymbolf\nF =\nbeginbmatrix\nf_1(x_1)  f_2(x_1)  cdots  f_N_d(x_1) \nf_1(x_2)  f_2(x_2)  cdots  f_N_d(x_2) \nvdots  vdots  ddots  vdots \nf_1(x_n)  f_2(x_n)  cdots  f_N_d(x_n)\nendbmatrixand polynomial vector f = beginbmatrix f_1(x_0)  f_2(x_0)  cdots  f_N_d(x_0) endbmatrix^top.The variogram matrix is constructed instead of the covariance matrix:newcommandGboldsymbolGamma\nnewcommandgboldsymbolgamma\nG =\nbeginbmatrix\ngamma(x_1x_1)  gamma(x_1x_2)  cdots  gamma(x_1x_n) \ngamma(x_2x_1)  gamma(x_2x_2)  cdots  gamma(x_2x_n) \nvdots  vdots  ddots  vdots \ngamma(x_nx_1)  gamma(x_nx_2)  cdots  gamma(x_nx_n)\nendbmatrixwith gamma(x_ix_j) = cov(x_0x_0) - cov(x_ix_j). The variogram is also evaluated at the estimation location g = beginbmatrix gamma(x_1x_0)  gamma(x_2x_0)  cdots  gamma(x_nx_0) endbmatrix^top.The resulting linear system is:beginbmatrix\nG  F \nF^top  boldsymbol0\nendbmatrix\nbeginbmatrix\nl \nboldsymbolnu\nendbmatrix\n=\nbeginbmatrix\ng \nf\nendbmatrixwith boldsymbolnu the Lagrange multipliers associated with the universal constraints. The mean and variance at location x_0 are given by:mu(x_0) = z^top lsigma^2(x_0) = beginbmatrixg  fendbmatrix^top beginbmatrixl  boldsymbolnuendbmatrixUniversalKriging"
},

{
    "location": "library.html#",
    "page": "Function Reference",
    "title": "Function Reference",
    "category": "page",
    "text": ""
},

{
    "location": "library.html#Types-1",
    "page": "Function Reference",
    "title": "Types",
    "category": "section",
    "text": "Modules = [GeoStats]\nOrder = [:type]"
},

{
    "location": "library.html#Functions-1",
    "page": "Function Reference",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [GeoStats]\nOrder = [:function]"
},

]}
