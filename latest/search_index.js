var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "(Image: GeoStatsLogo)(Image: Build Status) (Image: GeoStats) (Image: Coverage Status) (Image: Stable Documentation) (Image: Latest Documentation)"
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "High-performance implementations of geostatistical algorithms for the Julia programming language. This package is in its initial development, and currently only contains Kriging estimation methods. More features will be added as the Julia type system matures."
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
    "text": "Below is a quick example of usage:using GeoStats\nsrand(2017) # hide\n\n# create some data\ndim, nobs = 3, 10\nX = rand(dim, nobs)\nz = rand(nobs)\n\n# target location\nxₒ = rand(dim)\n\n# define a variogram model\nγ = GaussianVariogram(sill=1., range=1., nugget=0.)\n\n# define an estimator (i.e. build the Kriging system)\nsimkrig = SimpleKriging(X, z, γ, mean(z))\nordkrig = OrdinaryKriging(X, z, γ)\nunikrig = UniversalKriging(X, z, γ, 0)\n\n# estimate at target location\nμ, σ² = estimate(simkrig, xₒ)\nprintln(\"Simple Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide\nμ, σ² = estimate(ordkrig, xₒ)\nprintln(\"Ordinary Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide\nμ, σ² = estimate(unikrig, xₒ)\nprintln(\"Universal Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide"
},

{
    "location": "empirical_variograms.html#",
    "page": "Empirical variograms",
    "title": "Empirical variograms",
    "category": "page",
    "text": ""
},

{
    "location": "empirical_variograms.html#Empirical-variograms-1",
    "page": "Empirical variograms",
    "title": "Empirical variograms",
    "category": "section",
    "text": "An empirical variogram has the form:newcommandxboldsymbolx\nhatgamma(h) = frac12N(h) sum_(ij) in N(h) (z_i - z_j)^2where N(h) = left(ij) mid x_i - x_j = hright is the set of pairs of locations at a distance h and N(h) is the cardinality of the set.This package currently implements a simple ominidirectional variogram. Other more flexible types are planned for future releases.Empirical variograms can be estimated using general distance functions, please see Distance functions."
},

{
    "location": "empirical_variograms.html#GeoStats.EmpiricalVariogram",
    "page": "Empirical variograms",
    "title": "GeoStats.EmpiricalVariogram",
    "category": "Type",
    "text": "EmpiricalVariogram(X, z, [optional parameters])\n\nComputes the empirical (a.k.a. experimental) omnidirectional (semi-)variogram from data locations X and values z.\n\nEmpiricalVariogram(spatialdata, var, [optional parameters])\n\nAlternatively, compute the variogram for the variable var stored in a spatialdata object.\n\nParameters\n\nnbins - number of bins (default to 20)\nmaxlag - maximum lag distance (default to maximum lag of data)\ndistance - custom distance function\n\n\n\n"
},

{
    "location": "empirical_variograms.html#Omnidirectional-1",
    "page": "Empirical variograms",
    "title": "Omnidirectional",
    "category": "section",
    "text": "EmpiricalVariogram"
},

{
    "location": "theoretical_variograms.html#",
    "page": "Theoretical variograms",
    "title": "Theoretical variograms",
    "category": "page",
    "text": ""
},

{
    "location": "theoretical_variograms.html#GeoStats.isstationary",
    "page": "Theoretical variograms",
    "title": "GeoStats.isstationary",
    "category": "Function",
    "text": "isstationary(γ)\n\nCheck if variogram γ possesses the 2nd-order stationary property.\n\n\n\n"
},

{
    "location": "theoretical_variograms.html#Theoretical-variograms-1",
    "page": "Theoretical variograms",
    "title": "Theoretical variograms",
    "category": "section",
    "text": "newcommandxboldsymbolx\nnewcommandRmathbbR\nnewcommand1mathbb1In an intrinsic isotropic model, the variogram is only a function of the distance between any two points x_1x_2 in R^m:gamma(x_1x_2) = gamma(x_1 - x_2) = gamma(h)Under the additional assumption of 2nd-order stationarity, the well-known covariance is directly related via gamma(h) = cov(0) - cov(h). Anisotropic models are easily obtained by defining an ellipsoid distance in place of the Euclidean distance. For a list of available distances, please see Distance functions.This package implements a few commonly used and other more excentric variogram models. They all share the same default parameters:sill=1\nrange=1\nnugget=0\ndistance=EuclideanDistance()Some of them have extra parameters that can be set with keyword arguments:GaussianVariogram(nugget=.1) # set nugget effect\nMaternVariogram(order=1) # set order of Bessel functionAdditionally, a composite (additive) variogram model gamma(h) = gamma_1(h) + gamma_2(h) + cdots gamma_n(h) can be constructed from a list of variogram models:CompositeVariogram(GaussianVariogram(), ExponentialVariogram())Like the other variogram models, a composite variogram gamma can be evaluated as an isotropic model gamma(h) or as a model with a custom distance implicitly defined by taking into account its individual components gamma(x_1x_2).Finally, the 2nd-order stationarity property of a variogram can be checked with the isstationary method:isstationary"
},

{
    "location": "theoretical_variograms.html#GeoStats.GaussianVariogram",
    "page": "Theoretical variograms",
    "title": "GeoStats.GaussianVariogram",
    "category": "Type",
    "text": "GaussianVariogram(sill=s, range=r, nugget=n, distance=d)\n\nA Gaussian variogram with sill s, range r and nugget n. Optionally, use a custom distance d.\n\n\n\n"
},

{
    "location": "theoretical_variograms.html#Gaussian-1",
    "page": "Theoretical variograms",
    "title": "Gaussian",
    "category": "section",
    "text": "gamma(h) = (s - n) left1 - expleft(-left(frachrright)^2right)right + n cdot 1_(0infty)(h)GaussianVariogram"
},

{
    "location": "theoretical_variograms.html#GeoStats.ExponentialVariogram",
    "page": "Theoretical variograms",
    "title": "GeoStats.ExponentialVariogram",
    "category": "Type",
    "text": "ExponentialVariogram(sill=s, range=r, nugget=n, distance=d)\n\nAn exponential variogram with sill s, range r and nugget n. Optionally, use a custom distance d.\n\n\n\n"
},

{
    "location": "theoretical_variograms.html#Exponential-1",
    "page": "Theoretical variograms",
    "title": "Exponential",
    "category": "section",
    "text": "gamma(h) = (s - n) left1 - expleft(-frachrright)right + n cdot 1_(0infty)(h)\nExponentialVariogram"
},

{
    "location": "theoretical_variograms.html#GeoStats.MaternVariogram",
    "page": "Theoretical variograms",
    "title": "GeoStats.MaternVariogram",
    "category": "Type",
    "text": "MaternVariogram(sill=s, range=r, nugget=n, order=ν, distance=d)\n\nA Matérn variogram with sill s, range r and nugget n. The parameter ν is the order of the Bessel function. Optionally, use a custom distance d.\n\n\n\n"
},

{
    "location": "theoretical_variograms.html#Matern-1",
    "page": "Theoretical variograms",
    "title": "Matern",
    "category": "section",
    "text": "gamma(h) = (s - n) left1 - frac2^1-nuGamma(nu) left(sqrt2nufrachrright)^nu K_nuleft(sqrt2nufrachrright)right + n cdot 1_(0infty)(h)MaternVariogram"
},

{
    "location": "theoretical_variograms.html#GeoStats.SphericalVariogram",
    "page": "Theoretical variograms",
    "title": "GeoStats.SphericalVariogram",
    "category": "Type",
    "text": "SphericalVariogram(sill=s, range=r, nugget=n, distance=d)\n\nA spherical variogram with sill s, range r and nugget n. Optionally, use a custom distance d.\n\n\n\n"
},

{
    "location": "theoretical_variograms.html#Spherical-1",
    "page": "Theoretical variograms",
    "title": "Spherical",
    "category": "section",
    "text": "gamma(h) = (s - n) leftleft(frac32left(frachrright) + frac12left(frachrright)^3right) cdot 1_(0r)(h) + 1_rinfty)(h)right + n cdot 1_(0infty)(h)SphericalVariogram"
},

{
    "location": "theoretical_variograms.html#GeoStats.CubicVariogram",
    "page": "Theoretical variograms",
    "title": "GeoStats.CubicVariogram",
    "category": "Type",
    "text": "CubicVariogram(sill=s, range=r, nugget=n, distance=d)\n\nA cubic variogram with sill s, range r and nugget n. Optionally, use a custom distance d.\n\n\n\n"
},

{
    "location": "theoretical_variograms.html#Cubic-1",
    "page": "Theoretical variograms",
    "title": "Cubic",
    "category": "section",
    "text": "gamma(h) = (s - n) leftleft(7left(frachrright)^2 - frac354left(frachrright)^3 + frac72left(frachrright)^5 - frac34left(frachrright)^7right) cdot 1_(0r)(h) + 1_rinfty)(h)right + n cdot 1_(0infty)(h)CubicVariogram"
},

{
    "location": "theoretical_variograms.html#GeoStats.PentasphericalVariogram",
    "page": "Theoretical variograms",
    "title": "GeoStats.PentasphericalVariogram",
    "category": "Type",
    "text": "PentasphericalVariogram\n\nA pentaspherical variogram with sill s, range r and nugget n. Optionally, use a custom distance d.\n\n\n\n"
},

{
    "location": "theoretical_variograms.html#Pentaspherical-1",
    "page": "Theoretical variograms",
    "title": "Pentaspherical",
    "category": "section",
    "text": "gamma(h) = (s - n) leftleft(frac158left(frachrright) - frac54left(frachrright)^3 + frac38left(frachrright)^5right) cdot 1_(0r)(h) + 1_rinfty)(h)right + n cdot 1_(0infty)(h)PentasphericalVariogram"
},

{
    "location": "theoretical_variograms.html#GeoStats.PowerVariogram",
    "page": "Theoretical variograms",
    "title": "GeoStats.PowerVariogram",
    "category": "Type",
    "text": "PowerVariogram(scaling=s, exponent=a, nugget=n, distance=d)\n\nA power variogram with scaling s, exponent a and nugget n. Optionally, use a custom distance d.\n\n\n\n"
},

{
    "location": "theoretical_variograms.html#Power-1",
    "page": "Theoretical variograms",
    "title": "Power",
    "category": "section",
    "text": "gamma(h) = sh^a + n cdot 1_(0infty)(h)PowerVariogram"
},

{
    "location": "theoretical_variograms.html#GeoStats.SineHoleVariogram",
    "page": "Theoretical variograms",
    "title": "GeoStats.SineHoleVariogram",
    "category": "Type",
    "text": "SineHoleVariogram(sill=s, range=r, nugget=n, distance=d)\n\nA sine hole variogram with sill s, range r and nugget n. Optionally, use a custom distance d.\n\n\n\n"
},

{
    "location": "theoretical_variograms.html#Sine-hole-1",
    "page": "Theoretical variograms",
    "title": "Sine hole",
    "category": "section",
    "text": "gamma(h) = (s - n) left1 - fracsin(pi h  r)pi h  rright + n cdot 1_(0infty)(h)SineHoleVariogram"
},

{
    "location": "theoretical_variograms.html#GeoStats.CompositeVariogram",
    "page": "Theoretical variograms",
    "title": "GeoStats.CompositeVariogram",
    "category": "Type",
    "text": "CompositeVariogram(γ₁, γ₂, ..., γₙ)\n\nA composite (additive) model of variograms γ(h) = γ₁(h) + γ₂(h) + ⋯ + γₙ(h).\n\n\n\n"
},

{
    "location": "theoretical_variograms.html#Composite-1",
    "page": "Theoretical variograms",
    "title": "Composite",
    "category": "section",
    "text": "CompositeVariogram"
},

{
    "location": "estimators.html#",
    "page": "Kriging estimators",
    "title": "Kriging estimators",
    "category": "page",
    "text": ""
},

{
    "location": "estimators.html#GeoStats.estimate",
    "page": "Kriging estimators",
    "title": "GeoStats.estimate",
    "category": "Function",
    "text": "estimate(estimator, xₒ)\n\nCompute mean and variance for the estimator at coordinates xₒ.\n\n\n\n"
},

{
    "location": "estimators.html#GeoStats.fit!",
    "page": "Kriging estimators",
    "title": "GeoStats.fit!",
    "category": "Function",
    "text": "fit!(estimator, X, z)\n\nBuild Kriging system from locations X with values z and save factorization in estimator.\n\n\n\n"
},

{
    "location": "estimators.html#GeoStats.weights",
    "page": "Kriging estimators",
    "title": "GeoStats.weights",
    "category": "Function",
    "text": "weights(estimator, xₒ)\n\nCompute the weights λ (and Lagrange multipliers ν) for the estimator at coordinates xₒ.\n\n\n\n"
},

{
    "location": "estimators.html#Kriging-estimators-1",
    "page": "Kriging estimators",
    "title": "Kriging estimators",
    "category": "section",
    "text": "A Kriging estimator has the form:newcommandxboldsymbolx\nnewcommandRmathbbR\nhatZ(x_0) = lambda_1 Z(x_1) + lambda_2 Z(x_2) + cdots + lambda_n Z(x_n)quad x_i in R^m lambda_i in Rwith Zcolon R^m times Omega to R a random field.This package implements the following Kriging variants:Simple Kriging\nOrdinary Kriging\nUniversal Kriging\nExternal Drift KrigingAll these variants follow the same interface: an estimator object is first created with a given data configuration and variogram model, and then estimates are made at various locations.The object construction takes care of building the Kriging system and factorizing the LHS with an appropriate decomposition (e.g. Cholesky, LU). The estimate method performs the estimation at a given location:estimateA typical use of the interface is as follows:# build and factorize the system\nsimkrig = SimpleKriging(X, z, γ, mean(z))\n\n# estimate at various locations\nfor xₒ in locations\n  μ, σ² = estimate(simkrig, xₒ)\nendIn case the data configuration needs to be changed in a loop (e.g. sequential Gaussian simulation), one can keep all the parameters fixed and only update the factorization with the fit! method:fit!For advanced users, the Kriging weights and Lagrange multipliers at a given location can be accessed with the weights method. This method returns an AbstractWeights object containing a field λ for the weights and a field ν for the Lagrange multipliers:weights# weights and Lagrange multipliers\nOKweights = weights(ordkrig, xₒ)\nOKweights.λ, OKweights.ν"
},

{
    "location": "estimators.html#GeoStats.SimpleKriging",
    "page": "Kriging estimators",
    "title": "GeoStats.SimpleKriging",
    "category": "Type",
    "text": "SimpleKriging(X, z, γ, μ)\n\nParameters\n\nX ∈ ℜ^(mxn) - matrix of data locations\nz ∈ ℜⁿ      - vector of observations for X\nγ           - variogram model\nμ ∈ ℜ       - mean of z\n\nNotes\n\nSimple Kriging requires stationary variograms\n\n\n\n"
},

{
    "location": "estimators.html#Simple-Kriging-1",
    "page": "Kriging estimators",
    "title": "Simple Kriging",
    "category": "section",
    "text": "In Simple Kriging, the mean mu of the random field is assumed to be constant and known. The resulting linear system is:newcommandCboldsymbolC\nnewcommandcboldsymbolc\nnewcommandlboldsymbollambda\nnewcommand1boldsymbol1\nnewcommandzboldsymbolz\nbeginbmatrix\ncov(x_1x_2)  cov(x_1x_2)  cdots  cov(x_1x_n) \ncov(x_2x_1)  cov(x_2x_2)  cdots  cov(x_2x_n) \nvdots  vdots  ddots  vdots \ncov(x_nx_1)  cov(x_nx_2)  cdots  cov(x_nx_n)\nendbmatrix\nbeginbmatrix\nlambda_1 \nlambda_2 \nvdots \nlambda_n\nendbmatrix\n=\nbeginbmatrix\ncov(x_1x_0) \ncov(x_2x_0) \nvdots \ncov(x_nx_0)\nendbmatrixor in matricial form Cl = c. We subtract the given mean from the observations boldsymboly = z - mu 1 and compute the mean and variance at location x_0:mu(x_0) = mu + boldsymboly^top lsigma^2(x_0) = cov(0) - c^top lSimpleKriging"
},

{
    "location": "estimators.html#GeoStats.OrdinaryKriging",
    "page": "Kriging estimators",
    "title": "GeoStats.OrdinaryKriging",
    "category": "Type",
    "text": "OrdinaryKriging(X, z, γ)\n\nParameters\n\nX ∈ ℜ^(mxn) - matrix of data locations\nz ∈ ℜⁿ      - vector of observations for X\nγ           - variogram model\n\n\n\n"
},

{
    "location": "estimators.html#Ordinary-Kriging-1",
    "page": "Kriging estimators",
    "title": "Ordinary Kriging",
    "category": "section",
    "text": "In Ordinary Kriging the mean of the random field is assumed to be constant and unknown. The resulting linear system is:newcommandGboldsymbolGamma\nnewcommandgboldsymbolgamma\nbeginbmatrix\nG  1 \n1^top  0\nendbmatrix\nbeginbmatrix\nl \nnu\nendbmatrix\n=\nbeginbmatrix\ng \n1\nendbmatrixwith nu the Lagrange multiplier associated with the constraint 1^top l = 1. The mean and variance at location x_0 are given by:mu(x_0) = z^top lambdasigma^2(x_0) =  beginbmatrix g  1 endbmatrix^top beginbmatrix l  nu endbmatrixOrdinaryKriging"
},

{
    "location": "estimators.html#GeoStats.UniversalKriging",
    "page": "Kriging estimators",
    "title": "GeoStats.UniversalKriging",
    "category": "Type",
    "text": "UniversalKriging(X, z, γ, degree)\n\nParameters\n\nX ∈ ℜ^(mxn) - matrix of data locations\nz ∈ ℜⁿ      - vector of observations for X\nγ           - variogram model\ndegree      - polynomial degree for the mean\n\nNotes\n\nOrdinaryKriging is recovered for 0th degree polynomial\nFor non-polynomial mean, see ExternalDriftKriging\n\n\n\n"
},

{
    "location": "estimators.html#Universal-Kriging-1",
    "page": "Kriging estimators",
    "title": "Universal Kriging",
    "category": "section",
    "text": "In Universal Kriging, the mean of the random field is assumed to be a polynomial of the spatial coordinates:mu(x) = sum_k=1^N_d beta_k f_k(x)with N_d monomials f_k of degree up to d. For example, in 2D there are 6 monomials of degree up to 2:mu(x_1x_2) =  beta_1 1 + beta_2 x_1 + beta_3 x_2 + beta_4 x_1 x_2 + beta_5 x_1^2 + beta_6 x_2^2The choice of the degree d determines the size of the polynomial matrixnewcommandFboldsymbolF\nnewcommandfboldsymbolf\nF =\nbeginbmatrix\nf_1(x_1)  f_2(x_1)  cdots  f_N_d(x_1) \nf_1(x_2)  f_2(x_2)  cdots  f_N_d(x_2) \nvdots  vdots  ddots  vdots \nf_1(x_n)  f_2(x_n)  cdots  f_N_d(x_n)\nendbmatrixand polynomial vector f = beginbmatrix f_1(x_0)  f_2(x_0)  cdots  f_N_d(x_0) endbmatrix^top.The variogram determines the variogram matrix:G =\nbeginbmatrix\ngamma(x_1x_1)  gamma(x_1x_2)  cdots  gamma(x_1x_n) \ngamma(x_2x_1)  gamma(x_2x_2)  cdots  gamma(x_2x_n) \nvdots  vdots  ddots  vdots \ngamma(x_nx_1)  gamma(x_nx_2)  cdots  gamma(x_nx_n)\nendbmatrixand the variogram vector g = beginbmatrix gamma(x_1x_0)  gamma(x_2x_0)  cdots  gamma(x_nx_0) endbmatrix^top.The resulting linear system is:beginbmatrix\nG  F \nF^top  boldsymbol0\nendbmatrix\nbeginbmatrix\nl \nboldsymbolnu\nendbmatrix\n=\nbeginbmatrix\ng \nf\nendbmatrixwith boldsymbolnu the Lagrange multipliers associated with the universal constraints. The mean and variance at location x_0 are given by:mu(x_0) = z^top lsigma^2(x_0) = beginbmatrixg  fendbmatrix^top beginbmatrixl  boldsymbolnuendbmatrixUniversalKriging"
},

{
    "location": "estimators.html#GeoStats.ExternalDriftKriging",
    "page": "Kriging estimators",
    "title": "GeoStats.ExternalDriftKriging",
    "category": "Type",
    "text": "ExternalDriftKriging(X, z, γ, drifts)\n\nParameters\n\nX ∈ ℜ^(mxn) - matrix of data locations\nz ∈ ℜⁿ      - vector of observations for X\nγ           - variogram model\ndrifts      - vector of external drift functions m: ℜᵐ ↦ ℜ\n\nNotes\n\nExternal drift functions should be smooth\nKriging system with external drift is often unstable\nInclude a constant drift (e.g. x->1) for unbiased estimation\nOrdinaryKriging is recovered for drifts = [x->1]\nFor polynomial mean, see UniversalKriging\n\n\n\n"
},

{
    "location": "estimators.html#External-Drift-Kriging-1",
    "page": "Kriging estimators",
    "title": "External Drift Kriging",
    "category": "section",
    "text": "In External Drift Kriging, the mean of the random field is assumed to be a combination of known smooth functions:mu(x) = sum_k beta_k m_k(x)Differently than Universal Kriging, the functions m_k are not necessarily polynomials of the spatial coordinates. In practice, they represent a list of variables that is strongly correlated (and co-located) with the variable being estimated.External drifts are known to cause numerical instability. Give preference to other Kriging variants if possible.ExternalDriftKriging"
},

{
    "location": "distances.html#",
    "page": "Distance functions",
    "title": "Distance functions",
    "category": "page",
    "text": ""
},

{
    "location": "distances.html#Distance-functions-1",
    "page": "Distance functions",
    "title": "Distance functions",
    "category": "section",
    "text": "A set of commonly used distance functions is provided in this package for use in geostatistical algorithms. They can be passed to variograms in order to:Model anisotropy (e.g. ellipsoid distance)\nPerform geostatistical simulation on non-Euclidean coordinate systems (e.g. haversine distance)Custom distance functions are particularly useful if 3D locations are projected on a 2D map by means of a non-trivial transformation. In this case, a geodesic distance can be defined to properly account for spatial distortions at large scales."
},

{
    "location": "distances.html#GeoStats.EuclideanDistance",
    "page": "Distance functions",
    "title": "GeoStats.EuclideanDistance",
    "category": "Type",
    "text": "EuclideanDistance\n\nThe Euclidean distance ||x-y||₂\n\n\n\n"
},

{
    "location": "distances.html#Euclidean-1",
    "page": "Distance functions",
    "title": "Euclidean",
    "category": "section",
    "text": "newcommandxboldsymbolx\nnewcommandyboldsymboly\nd(xy) = sqrt(x-y)^top (x-y)EuclideanDistance"
},

{
    "location": "distances.html#GeoStats.EllipsoidDistance",
    "page": "Distance functions",
    "title": "GeoStats.EllipsoidDistance",
    "category": "Type",
    "text": "EllipsoidDistance(semiaxes, angles)\n\nA distance defined by an ellipsoid with given semiaxes and rotation angles.\n\nFor 2D ellipsoids, there are two semiaxes and one rotation angle.\nFor 3D ellipsoids, there are three semiaxes and three rotation angles.\n\nExamples\n\n2D ellipsoid making 45ᵒ with the horizontal axis:\n\njulia> EllipsoidDistance([1.0,0.5], [π/2])\n\n3D ellipsoid rotated by 45ᵒ in the xy plane:\n\njulia> EllipsoidDistance([1.0,0.5,0.5], [π/2,0.0,0.0])\n\nNotes\n\nThe positive definite matrix representing the ellipsoid is assembled once during object construction and cached for fast evaluation.\n\n\n\n"
},

{
    "location": "distances.html#Ellipsoid-1",
    "page": "Distance functions",
    "title": "Ellipsoid",
    "category": "section",
    "text": "The ellipsoid distance can be used to model anisotropy. The semiaxes of the ellipsoid represent correlation lengths that can be rotated and aligned with target directions.d(xy) = sqrt(x-y)^top boldsymbolA (x-y)EllipsoidDistance"
},

{
    "location": "distances.html#GeoStats.HaversineDistance",
    "page": "Distance functions",
    "title": "GeoStats.HaversineDistance",
    "category": "Type",
    "text": "HaversineDistance(radius)\n\nThe haversine distance between two locations on a sphere of given radius.\n\nLocations are described with longitude and latitude in degrees and the radius of the Earth is used by default (≈ 6371km). The computed distance has the same units as that of the radius.\n\nNotes\n\nThe haversine formula is widely used to approximate the geodesic distance between two points at the surface of the Earth. The error from approximating the Earth as a sphere is typically negligible for most applications. It is no more than 0.3%.\n\n\n\n"
},

{
    "location": "distances.html#Haversine-1",
    "page": "Distance functions",
    "title": "Haversine",
    "category": "section",
    "text": "The haversine distance can be used to perform geostatistical simulation directly on a sphere. It approximates the geodesic distance between two pairs of latitude/longitude.HaversineDistance"
},

{
    "location": "examples.html#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples.html#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "A set of Jupyter notebooks demonstrating the current functionality of the package is available in the examples folder. These notebooks are distributed with GeoStats.jl and can be run locally with GeoStats.examples().Want to contribute an example? Please check Contributing before submitting a pull request."
},

{
    "location": "plotting.html#",
    "page": "Plotting",
    "title": "Plotting",
    "category": "page",
    "text": ""
},

{
    "location": "plotting.html#Plotting-1",
    "page": "Plotting",
    "title": "Plotting",
    "category": "section",
    "text": "GeoStats.jl is integrated with the Julia Plots.jl API. This means that many objects defined in the package can be plotted directly without data format conversions.For example, below we plot various theoretical variograms with the plot command from Plots.jl:using GeoStats\nusing Plots\ngr(size=(600,400)) # hide\n\nplot(GaussianVariogram(), maxlag=3., label=\"Gaussian\")\nplot!(ExponentialVariogram(), maxlag=3., label=\"Exponential\")\nplot!(SphericalVariogram(), maxlag=3., label=\"Spherical\")\nplot!(MaternVariogram(), maxlag=3., label=\"Matern\")\npng(\"images/variograms.png\") # hide(Image: )"
},

{
    "location": "library.html#",
    "page": "Library",
    "title": "Library",
    "category": "page",
    "text": ""
},

{
    "location": "library.html#Types-1",
    "page": "Library",
    "title": "Types",
    "category": "section",
    "text": "Modules = [GeoStats]\nOrder = [:type]"
},

{
    "location": "library.html#Functions-1",
    "page": "Library",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [GeoStats]\nOrder = [:function]"
},

{
    "location": "contributing.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "contributing.html#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "First off, thank you for considering contributing to GeoStats.jl. It’s people like you that make this project so much fun."
},

{
    "location": "contributing.html#Reporting-issues-1",
    "page": "Contributing",
    "title": "Reporting issues",
    "category": "section",
    "text": "If you are experiencing issues or have discovered a bug, please report it on GitHub. To make the resolution process easier, please include the version of Julia and GeoStats.jl in your writeup. These can be found with two commands:julia> versioninfo()\njulia> Pkg.status(\"GeoStats\")"
},

{
    "location": "contributing.html#Feature-requests-1",
    "page": "Contributing",
    "title": "Feature requests",
    "category": "section",
    "text": "If you have suggestions of improvement or algorithms that you would like to see implemented in GeoStats.jl, please open an issue on GitHub. Suggestions as well as feature requests are very welcome."
},

{
    "location": "contributing.html#Code-contribution-1",
    "page": "Contributing",
    "title": "Code contribution",
    "category": "section",
    "text": "If you have code that you would like to contribute to GeoStats.jl, that is awesome! Please open an issue before you create the pull request on GitHub so that we make sure your idea is aligned with our goals for the project."
},

{
    "location": "about/author.html#",
    "page": "Author",
    "title": "Author",
    "category": "page",
    "text": "Júlio Hoffimann MendesI am a Ph.D. candidate in the Department of Energy Resources Engineering at Stanford University. You can find more about my research on my website. Below are some ways that we can connect:ResearchGate\nLinkedIn\nGitHub"
},

{
    "location": "about/license.html#",
    "page": "License",
    "title": "License",
    "category": "page",
    "text": "The GeoStats.jl package is licensed under the ISC License:Copyright (c) 2015, Júlio Hoffimann Mendes <juliohm@stanford.edu>\n\nPermission to use, copy, modify, and/or distribute this software for any\npurpose with or without fee is hereby granted, provided that the above\ncopyright notice and this permission notice appear in all copies.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\" AND THE AUTHOR DISCLAIMS ALL WARRANTIES\nWITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF\nMERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR\nANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES\nWHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN\nACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF\nOR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE."
},

]}
