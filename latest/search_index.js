var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "(Image: GeoStatsLogo)(Image: Build Status) (Image: GeoStats) (Image: Coverage Status) (Image: Stable Documentation) (Image: Latest Documentation)High-performance implementations of geostatistical algorithms for the Julia programming language."
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "Geostatistics (a.k.a. spatial statistics) is the branch of statistics that deals with spatial data. In many fields of science, such as mining engineering, hidrogeology, petroleum engineering, and environmental sciences, traditional regression techniques fail to capture spatiotemporal correlation, and therefore are not satisfactory tools for decision making.GeoStats.jl is an attempt to bring together bleeding-edge research in the geostatistics community into a comprehensive framework, and to empower researchers and practioners with a toolkit for fast assessment of different modeling approaches.The design of this package is the result of many years developing geostatistical software. I hope that it can serve to promote more collaboration between geostatisticians around the globe and to standardize this incredible science.If you would like to help support the project, please star the repository on GitHub and share it with your colleagues. If you are a developer, please check GeoStatsBase.jl."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Get the latest stable release with Julia's package manager:Pkg.add(\"GeoStats\")"
},

{
    "location": "index.html#Project-organization-1",
    "page": "Home",
    "title": "Project organization",
    "category": "section",
    "text": "The project is split into various packages:Package Description\nGeoStats.jl Main package containing Kriging-based solvers, and other geostatistical tools.\nGeoStatsImages.jl Training images for multiple-point geostatistical simulation.\nGslibIO.jl Utilities to read/write extended GSLIB files.\nGeoStatsBase.jl Base package containing problem and solution specifications (for developers).\nGeoStatsDevTools.jl Developer tools for writing new solvers (for developers).The main package (i.e. GeoStats.jl) is self-contained, and provides high-performance Kriging-based estimation/simulation algorithms over arbitrary domains. Other packages can be installed from the list above for additional functionality."
},

{
    "location": "index.html#Quick-example-1",
    "page": "Home",
    "title": "Quick example",
    "category": "section",
    "text": "Below is a quick preview of the high-level API, for the full example, please see Examples.using GeoStats\nusing Plots\n\n# data.csv:\n#    x,    y,       station, precipitation\n# 25.0, 25.0,     palo alto,           1.0\n# 50.0, 75.0,  redwood city,           0.0\n# 75.0, 50.0, mountain view,           1.0\n\n# read spreadsheet file containing spatial data\ngeodata = readtable(\"data.csv\", coordnames=[:x,:y])\n\n# define spatial domain (e.g. regular grid, point collection)\ngrid = RegularGrid{Float64}(100, 100)\n\n# define estimation problem for any data column(s) (e.g. :precipitation)\nproblem = EstimationProblem(geodata, grid, :precipitation)\n\n# solve the problem with any solver\nsolution = solve(problem, Kriging())\n\n# plot the solution\nplot(solution)(Image: EstimationSolution)"
},

{
    "location": "index.html#Low-level-API-1",
    "page": "Home",
    "title": "Low-level API",
    "category": "section",
    "text": "If you are interested in finer control, Kriging estimators can also be used directly:using GeoStats\nsrand(2017) # hide\n\n# create some data\ndim, nobs = 3, 10\nX = rand(dim, nobs)\nz = rand(nobs)\n\n# target location\nxₒ = rand(dim)\n\n# define a variogram model\nγ = GaussianVariogram(sill=1., range=1., nugget=0.)\n\n# define an estimator (i.e. build the Kriging system)\nsimkrig = SimpleKriging(X, z, γ, mean(z))\nordkrig = OrdinaryKriging(X, z, γ)\nunikrig = UniversalKriging(X, z, γ, 0)\n\n# estimate at target location\nμ, σ² = estimate(simkrig, xₒ)\nprintln(\"Simple Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide\nμ, σ² = estimate(ordkrig, xₒ)\nprintln(\"Ordinary Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide\nμ, σ² = estimate(unikrig, xₒ)\nprintln(\"Universal Kriging:\") # hide\nprintln(\"  μ = $μ, σ² = $σ²\") # hide"
},

{
    "location": "problems_and_solvers.html#",
    "page": "Problems and solvers",
    "title": "Problems and solvers",
    "category": "page",
    "text": ""
},

{
    "location": "problems_and_solvers.html#Problems-and-solvers-1",
    "page": "Problems and solvers",
    "title": "Problems and solvers",
    "category": "section",
    "text": "One of the greatest features of GeoStats.jl is the ability to define geostatistical problems independently of the solution strategy. This design allows researchers and practioners to perform fair comparisons between different solvers. It is perhaps the single most important contribution of this project.If you are an experienced user of geostatistics or if you do research in the field, you know how hard it is to compare algorithms fairly. Often a new algorithm is proposed in the literature, and yet the task of comparing it with the state of the art is quite demanding. Even when a comparison is made by the author after a great amount of effort, it is inevitably biased.Part of this issue is attributed to the fact that a general definition of the problem is missing. What is it that we call an \"estimation problem\" in geostatistics? What about \"stochastic simulation\"? The answer to these questions is given below in the form of code."
},

{
    "location": "problems_and_solvers.html#GeoStatsBase.EstimationProblem",
    "page": "Problems and solvers",
    "title": "GeoStatsBase.EstimationProblem",
    "category": "Type",
    "text": "EstimationProblem(spatialdata, domain, targetvars)\n\nA spatial estimation problem on a given domain in which the variables to be estimated are listed in targetvars. The data of the problem is stored in spatialdata.\n\nExamples\n\nCreate an estimation problem for rainfall precipitation measurements:\n\njulia> EstimationProblem(spatialdata, domain, :precipitation)\n\nCreate an estimation problem for precipitation and CO₂:\n\njulia> EstimationProblem(spatialdata, domain, [:precipitation, :CO₂])\n\n\n\n"
},

{
    "location": "problems_and_solvers.html#Estimation-1",
    "page": "Problems and solvers",
    "title": "Estimation",
    "category": "section",
    "text": "An estimation problem in geostatitsics is a triplet:Spatial data (i.e. data with coordinates)\nSpatial domain (e.g. regular grid, point collection)\nTarget variables (or variables to be estimated)Each of these components is constructed separately, and then grouped (no memory is copied) in an EstimationProblem.EstimationProblemPlease check Spatial data and Domains for currently implemented data and domain types."
},

{
    "location": "problems_and_solvers.html#GeoStatsBase.SimulationProblem",
    "page": "Problems and solvers",
    "title": "GeoStatsBase.SimulationProblem",
    "category": "Type",
    "text": "SimulationProblem(spatialdata, domain, targetvars, nreals)\nSimulationProblem(domain, targetvars, nreals)\n\nA spatial simulation problem on a given domain in which the variables to be simulated are listed in targetvars.\n\nFor conditional simulation, the data of the problem is stored in spatialdata.\n\nFor unconditional simulation, a dictionary targetvars must be provided mapping variable names to their types.\n\nIn both cases, a number nreals of realizations is requested.\n\nExamples\n\nCreate a conditional simulation problem for porosity and permeability with 100 realizations:\n\njulia> SimulationProblem(spatialdata, domain, [:porosity,:permeability], 100)\n\nCreate an unconditional simulation problem for porosity and facies type with 100 realizations:\n\njulia> SimulationProblem(domain, Dict(:porosity => Float64, :facies => Int), 100)\n\nNotes\n\nTo check if a simulation problem has data (i.e. conditional vs. unconditional) use the hasdata method.\n\n\n\n"
},

{
    "location": "problems_and_solvers.html#GeoStatsBase.hasdata",
    "page": "Problems and solvers",
    "title": "GeoStatsBase.hasdata",
    "category": "Function",
    "text": "hasdata(problem)\n\nReturn true if simulation problem has data.\n\n\n\n"
},

{
    "location": "problems_and_solvers.html#Simulation-1",
    "page": "Problems and solvers",
    "title": "Simulation",
    "category": "section",
    "text": "Likewise, a stochastic simulation problem in geostatistics is represented with the same triplet. However, the spatial data in this case is optional in order to accomodate the concept of conditional versus unconditional simulation.SimulationProblemhasdata"
},

{
    "location": "problems_and_solvers.html#Solvers-1",
    "page": "Problems and solvers",
    "title": "Solvers",
    "category": "section",
    "text": "Below is the list of solvers distributed with GeoStats.jl. For more solvers, please check the main project page."
},

{
    "location": "problems_and_solvers.html#GeoStats.Kriging",
    "page": "Problems and solvers",
    "title": "GeoStats.Kriging",
    "category": "Type",
    "text": "Kriging(var₁=>param₁, var₂=>param₂, ...)\n\nA polyalgorithm Kriging estimation solver.\n\nEach pair var=>param specifies the KrigingParam param for the Kriging variable var. In order to avoid boilerplate code, the constructor expects pairs of Symbol and NamedTuple instead.\n\nParameters\n\nvariogram - Variogram model (default to GaussianVariogram())\nmean      - Simple Kriging mean\ndegree    - Universal Kriging degree\ndrifts    - External Drift Kriging drift functions\n\nLatter options override former options. For example, by specifying drifts, the user is telling the algorithm to ignore degree and mean. If no option is specified, Ordinary Kriging is used by default with the variogram only.\n\nExamples\n\nSolve the variable :var₁ with Simple Kriging by specifying the mean, and the variable :var₂ with Universal Kriging by specifying the degree and the variogram model.\n\njulia> Kriging(\n  :var₁ => @NT(mean=1.),\n  :var₂ => @NT(degree=1, variogram=SphericalVariogram(range=20.))\n)\n\nSolve all variables of the problem with the default parameters (i.e. Ordinary Kriging with unit Gaussian variogram):\n\njulia> Kriging()\n\nNotes\n\nThe prefix @NT extends for NamedTuple. It won't be necessary in Julia v0.7 and beyond.\n\n\n\n\n\n"
},

{
    "location": "problems_and_solvers.html#Estimation-2",
    "page": "Problems and solvers",
    "title": "Estimation",
    "category": "section",
    "text": "Kriging"
},

{
    "location": "problems_and_solvers.html#GeoStats.SeqGaussSim",
    "page": "Problems and solvers",
    "title": "GeoStats.SeqGaussSim",
    "category": "Type",
    "text": "SeqGaussSim(var₁=>param₁, var₂=>param₂, ...)\n\nA polyalgorithm sequential Gaussian simulation solver.\n\nEach pair var=>param specifies the SeqGaussSimParam param for the simulation variable var. In order to avoid boilerplate code, the constructor expects pairs of Symbol and NamedTuple instead. See Kriging documentation for examples.\n\nParameters\n\nvariogram - Variogram model (default to GaussianVariogram())\nmean      - Simple Kriging mean\ndegree    - Universal Kriging degree\ndrifts    - External Drift Kriging drift functions\n\nLatter options override former options. For example, by specifying drifts, the user is telling the algorithm to ignore degree and mean. If no option is specified, Ordinary Kriging is used by default with the variogram only.\n\npath         - Simulation path (default to :random)\nneighradius  - Radius of search neighborhood (default to 10.)\nmaxneighbors - Maximum number of neighbors (default to 10)\n\n\n\n\n\n"
},

{
    "location": "problems_and_solvers.html#Simulation-2",
    "page": "Problems and solvers",
    "title": "Simulation",
    "category": "section",
    "text": "SeqGaussSim"
},

{
    "location": "spatialdata.html#",
    "page": "Spatial data",
    "title": "Spatial data",
    "category": "page",
    "text": ""
},

{
    "location": "spatialdata.html#Spatial-data-1",
    "page": "Spatial data",
    "title": "Spatial data",
    "category": "section",
    "text": "In GeoStats.jl, data and domain are disconnected one from another. This design enables levels of parallelism that would be difficult or even impossible to implement otherwise.The data is accessed by solvers only if strictly necessary. One of our goals is to be able to handle massive datasets that may not fit in random-access memory (RAM)."
},

{
    "location": "spatialdata.html#GeoStatsDevTools.GeoDataFrame",
    "page": "Spatial data",
    "title": "GeoStatsDevTools.GeoDataFrame",
    "category": "Type",
    "text": "GeoDataFrame(data, coordnames)\n\nA dataframe object data with additional metadata for tracking the columns coordnames that represent spatial coordinates.\n\nExamples\n\nIf the data was already loaded in a normal DataFrame data, and there exists columns named x, y and z, wrap the data and specify the column names:\n\njulia> GeoDataFrame(data, [:x,:y,:z])\n\nAlternatively, load the data directly into a GeoDataFrame object by using the method readtable.\n\nNotes\n\nThis type is a lightweight wrapper over Julia's DataFrame types. No additional storage is required other than a vector of symbols with the columns names representing spatial coordinates.\n\n\n\n"
},

{
    "location": "spatialdata.html#GeoStatsDevTools.readtable",
    "page": "Spatial data",
    "title": "GeoStatsDevTools.readtable",
    "category": "Function",
    "text": "readtable(args; coordnames=[:x,:y,:z], kwargs)\n\nRead data from disk using DataFrames.readtable, optionally specifying the columns coordnames with spatial coordinates.\n\nThe arguments args and keyword arguments kwargs are forwarded to the DataFrames.readtable function, please check their documentation for more details.\n\nThis function returns a GeoDataFrame object.\n\n\n\n"
},

{
    "location": "spatialdata.html#GeoDataFrame-1",
    "page": "Spatial data",
    "title": "GeoDataFrame",
    "category": "section",
    "text": "For point (or hard) data in spreadsheet format (e.g. CSV, TSV), the GeoDataFrame object is a lightweight wrapper over Julia's DataFrame types.GeoDataFramereadtable"
},

{
    "location": "domains.html#",
    "page": "Domains",
    "title": "Domains",
    "category": "page",
    "text": ""
},

{
    "location": "domains.html#Domains-1",
    "page": "Domains",
    "title": "Domains",
    "category": "section",
    "text": "Although estimators such as Kriging do not require equally spaced samples, many software packages for geostatistical estimation/simulation restrict their implementations to regular grids. This is a big limitation, however; given that sampling campaigns and resource exploration is rarely regular.In GeoStats.jl, estimation/simulation can be performed on arbitrary domain types such as simple point collections, unstructured meshes, and tesselations. These types are implemented \"analytically\" in order to minimize memory access and maximize performance. In a regular grid, for example, the coordinates of a given location are known, and one can perform search without ever constructing a grid.Below is the list of currently implemented domain types. More options will be available in future releases."
},

{
    "location": "domains.html#GeoStatsDevTools.RegularGrid",
    "page": "Domains",
    "title": "GeoStatsDevTools.RegularGrid",
    "category": "Type",
    "text": "RegularGrid(dims, origin, spacing)\nRegularGrid{T}(dims)\n\nA regular grid with dimensions dims, lower left corner at origin and cell spacing spacing. The three arguments must have the same length.\n\nIn the first constructor, all the arguments are specified as vectors. In the second constructor, one needs to specify the type of the coordinates and the dimensions of the grid. In that case, the origin and spacing default to (0,0,...) and (1,1,...), respectively.\n\nExamples\n\nCreate a 3D regular grid with 100x100x50 locations:\n\njulia> RegularGrid{Float64}(100,100,50)\n\nCreate a 2D grid with 100x100 locations and origin at (10.,20.) units:\n\njulia> RegularGrid([100,100],[10.,20.],[1.,1.])\n\nNotes\n\nInternally, the vectors that are passed as arguments are converted into tuples that are stored in the stack instead of in the memory heap.\n\n\n\n"
},

{
    "location": "domains.html#Regular-grid-1",
    "page": "Domains",
    "title": "Regular grid",
    "category": "section",
    "text": "RegularGrid"
},

{
    "location": "domains.html#GeoStatsDevTools.PointCollection",
    "page": "Domains",
    "title": "GeoStatsDevTools.PointCollection",
    "category": "Type",
    "text": "PointCollection(coords)\n\nA collection of points with coordinate matrix coords. The number of rows is the dimensionality of the domain whereas the number of columns is the number of points.\n\n\n\n"
},

{
    "location": "domains.html#Point-collection-1",
    "page": "Domains",
    "title": "Point collection",
    "category": "section",
    "text": "PointCollection"
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
    "text": "A set of Jupyter notebooks demonstrating the current functionality of the project is available in the examples folder. These notebooks are distributed with GeoStats.jl and can be run locally with GeoStats.examples().Want to contribute an example? Please check the Contributing page before submitting a pull request. Contributions are very welcome, specially if they educate other users."
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
    "text": "First off, thank you for considering contributing to GeoStats.jl. It’s people like you that make this project so much fun. Below are a few suggestions to speed up the collaboration process."
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
