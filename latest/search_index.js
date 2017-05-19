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
    "text": "This package provides efficient geostatistical methods for the Julia programming language. It is in its initial development, and currently only implements Kriging estimation methods. More features will be added as the Julia type system matures."
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
    "text": "Below is a quick example of usage:using GeoStats\nsrand(2017) # hide\n\n# create some data\ndim, nobs = 3, 10\nX = rand(dim, nobs); z = rand(nobs)\n\n# target location\nxₒ = rand(dim)\n\n# define a covariance model\ncov = GaussianCovariance(0.,1.,1.) # nugget, sill and range\n\n# define an estimator (i.e. build the Kriging system)\nsimkrig = SimpleKriging(X, z, cov, mean(z))\nordkrig = OrdinaryKriging(X, z, cov)\nunikrig = UniversalKriging(X, z, cov, 1)\n\n# estimate at target location\nμ, σ² = estimate(simkrig, xₒ)\nprintln(\"Simple Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide\nμ, σ² = estimate(ordkrig, xₒ)\nprintln(\"Ordinary Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide\nμ, σ² = estimate(unikrig, xₒ)\nprintln(\"Universal Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide"
},

{
    "location": "estimation.html#",
    "page": "Estimation",
    "title": "Estimation",
    "category": "page",
    "text": ""
},

{
    "location": "estimation.html#Estimation-1",
    "page": "Estimation",
    "title": "Estimation",
    "category": "section",
    "text": "A Kriging estimator has the form:newcommandxboldsymbolx\nhatZ(x_0) = lambda_1 Z(x_1) + lambda_2 Z(x_2) + cdots + lambda_n Z(x_n)quad x_i in mathbbR^m lambda_i in mathbbRwith Zcolon mathbbR^m times Omega to mathbbR a random field.This package implements the following Kriging variants:Simple Kriging\nOrdinary Kriging\nUniversal Kriging (polynomial drift for the mean)All these variants follow the same interface: an estimator object is first created with a given data configuration and covariance model, and then estimates are made at various locations.The object construction takes care of building the Kriging system and factorizing the LHS with an appropriate decomposition (e.g. Cholesky, LU). The estimate method performs the estimation at a given location:# build and factorize the system\nsimkrig = SimpleKriging(X, z, cov, mean(z))\n\n# estimate at various locations\nfor xₒ in locations\n  μ, σ² = estimate(simkrig, xₒ)\nendIn case the data configuration needs to be changed in a loop (e.g. sequential Gaussian simulation), one can keep all the parameters fixed and only update the factorization with the fit! method:fit!(simkrig, Xnew, znew)"
},

{
    "location": "estimation.html#GeoStats.SimpleKriging",
    "page": "Estimation",
    "title": "GeoStats.SimpleKriging",
    "category": "Type",
    "text": "SimpleKriging(X, z, cov, μ)\n\nINPUTS:\n\n* X ∈ ℜ^(mxn) - matrix of data locations\n* z ∈ ℜⁿ      - vector of observations for X\n* cov         - covariance model\n* μ ∈ ℜ       - mean of z\n\n\n\n"
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
    "text": "OrdinaryKriging(X, z, cov)\n\nINPUTS:\n\n* X ∈ ℜ^(mxn) - matrix of data locations\n* z ∈ ℜⁿ      - vector of observations for X\n* cov         - covariance model\n\n\n\n"
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
    "text": "UniversalKriging(X, z, cov, degree)\n\nINPUTS:\n\n* X ∈ ℜ^(mxn) - matrix of data locations\n* z ∈ ℜⁿ      - vector of observations for X\n* cov         - covariance model\n* degree      - polynomial degree for the mean\n\nOrdinary Kriging is recovered for 0th degree polynomial.\n\n\n\n"
},

{
    "location": "estimation.html#Universal-Kriging-1",
    "page": "Estimation",
    "title": "Universal Kriging",
    "category": "section",
    "text": "In Universal Kriging, the mean of the random field is assumed to be a polynomial:mu(x) = sum_k=1^N_d beta_k f_k(x)with N_d monomials f_k of degree up to d. For example, in 2D there are 6 monomials of degree up to 2:mu(x_1x_2) =  beta_1 1 + beta_2 x_1 + beta_3 x_2 + beta_4 x_1 x_2 + beta_5 x_1^2 + beta_6 x_2^2The choice of the degree d determines the size of the polynomial matrixnewcommandFboldsymbolF\nnewcommandfboldsymbolf\nF =\nbeginbmatrix\nf_1(x_1)  f_2(x_1)  cdots  f_N_d(x_1) \nf_1(x_2)  f_2(x_2)  cdots  f_N_d(x_2) \nvdots  vdots  ddots  vdots \nf_1(x_n)  f_2(x_n)  cdots  f_N_d(x_n)\nendbmatrixand polynomial vector f = beginbmatrix f_1(x_0)  f_2(x_0)  cdots  f_N_d(x_0) endbmatrix^top.The variogram matrix is constructed instead of the covariance matrix:newcommandGboldsymbolGamma\nnewcommandgboldsymbolgamma\nG =\nbeginbmatrix\ngamma(x_1x_1)  gamma(x_1x_2)  cdots  gamma(x_1x_n) \ngamma(x_2x_1)  gamma(x_2x_2)  cdots  gamma(x_2x_n) \nvdots  vdots  ddots  vdots \ngamma(x_nx_1)  gamma(x_nx_2)  cdots  gamma(x_nx_n)\nendbmatrixwith gamma(x_ix_j) = cov(x_0x_0) - cov(x_ix_j). The variogram is also evaluated at the estimation location g = beginbmatrix gamma(x_1x_0)  gamma(x_2x_0)  cdots  gamma(x_nx_0) endbmatrix^top.The resulting linear system is:beginbmatrix\nG  F \nF^top  boldsymbol0\nendbmatrix\nbeginbmatrix\nl \nboldsymbolnu\nendbmatrix\n=\nbeginbmatrix\ng \nf\nendbmatrixwith boldsymbolnu the Lagrange multipliers associated with the universal constraints. The mean and variance at location x_0 are given by:mu(x_0) = z^top lsigma^2(x_0) = beginbmatrixg  fendbmatrix^top beginbmatrixl  boldsymbolnuendbmatrixUniversalKriging"
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
    "text": "newcommandxboldsymbolxIn a stationary isotropic model, the variogram is only a function of the distance between any two points x_1x_2 in mathbbR^m:gamma(x_1x_2) = gamma(x_1 - x_2) = gamma(h)The same holds for the covariance, which is directly related gamma(h) = cov(0) - cov(h). This package implements a few commonly used stationary models:"
},

{
    "location": "variograms.html#GeoStats.GaussianCovariance",
    "page": "Variograms",
    "title": "GeoStats.GaussianCovariance",
    "category": "Type",
    "text": "GaussianCovariance(n, s, r)\n\nINPUTS:\n\n* n ∈ ℜ - nugget\n* s ∈ ℜ - sill\n* r ∈ ℜ - range\n\n\n\n"
},

{
    "location": "variograms.html#Gaussian-1",
    "page": "Variograms",
    "title": "Gaussian",
    "category": "section",
    "text": "cov(h) = (s - n) cdot expleft(-left(frachrright)^2right)GaussianCovariance"
},

{
    "location": "variograms.html#GeoStats.SphericalCovariance",
    "page": "Variograms",
    "title": "GeoStats.SphericalCovariance",
    "category": "Type",
    "text": "SphericalCovariance(n, s, r)\n\nINPUTS:\n\n* n ∈ ℜ - nugget\n* s ∈ ℜ - sill\n* r ∈ ℜ - range\n\n\n\n"
},

{
    "location": "variograms.html#Spherical-1",
    "page": "Variograms",
    "title": "Spherical",
    "category": "section",
    "text": "cov(h) =\nbegincases\n(s - n) (1 - frac32left(frachrright) + frac12left(frachrright)^3)  textif  h leq r \n0  textotherwise\nendcasesSphericalCovariance"
},

{
    "location": "variograms.html#GeoStats.ExponentialCovariance",
    "page": "Variograms",
    "title": "GeoStats.ExponentialCovariance",
    "category": "Type",
    "text": "ExponentialCovariance(n, s, r)\n\nINPUTS:\n\n* n ∈ ℜ - nugget\n* s ∈ ℜ - sill\n* r ∈ ℜ - range\n\n\n\n"
},

{
    "location": "variograms.html#Exponential-1",
    "page": "Variograms",
    "title": "Exponential",
    "category": "section",
    "text": "cov(h) = (s - n) cdot expleft(-frachrright)ExponentialCovariance"
},

{
    "location": "library.html#",
    "page": "Function Reference",
    "title": "Function Reference",
    "category": "page",
    "text": ""
},

]}
