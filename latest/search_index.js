var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Overview-1",
    "page": "Home",
    "title": "Overview",
    "category": "section",
    "text": "Geostatistics in Julia.(Image: Build Status) (Image: GeoStats) (Image: Coverage Status)This project provides efficient geostatistical methods for the Julia programming language. It is in its initial development, and currently only implements Kriging estimation methods. More features will be added as the Julia type system matures."
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
    "text": "Below is a quick example of usage:using GeoStats\n\n# create some data\ndim, nobs = 3, 10\nX = rand(dim, nobs); z = rand(nobs)\n\n# target location\nxₒ = rand(dim)\n\n# define a covariance model\ncov = GaussianCovariance(0.,1.,1.) # nugget, sill and range\n\n# define an estimator (i.e. build the Kriging system)\nsimkrig = SimpleKriging(X, z, cov, mean(z))\nordkrig = OrdinaryKriging(X, z, cov)\nunikrig = UniversalKriging(X, z, cov, 1)\n\n# estimate at target location\nμ, σ² = estimate(simkrig, xₒ)\nμ, σ² = estimate(ordkrig, xₒ)\nμ, σ² = estimate(unikrig, xₒ)"
},

]}
