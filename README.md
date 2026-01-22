## Installation

First install the [`devtools`](https://rawgit.com/rstudio/cheatsheets/master/package-development.pdf) package and load it:

```R
install.packages("devtools")
library(devtools)
```

Then, you can install `clewsbiosphere` by directly loading it from github:

```R
devtools::install_github(repo = "https://github.com/stenzelf/clewsbiosphere")
library(clewsbiosphere)
```

or download from github and install from that folder
```bash
git clone git@github.com:stenzelf/clewsbiosphere.git /home/stenzel/clewsbiosphere
```

```R
devtools::load_all("/home/stenzel/clewsbiosphere/")
```

## Examples
The clews package provides functions:

`read_netcdf()`
`plot_lon_lat_array()`
`climate_plot()`
`read_in_zhang_file()`
`add_modelled_npp_to_zhang_npp_measurements()`
`scatter_plot()`

Base R functions you might need are:
```R
apply(array, keepTheseDimensions, functionToApply)
```
e.g. with apply(pre,c(1,2,4),sum) or apply(tmp,c(1,2,4),mean)

```R
pmin(vector1, vector2)
exp(x)
```
If you have limited RAM, you can remove variables that are not required anymore and clean the workspace
by these commands:
```R
rm(variable)
gc()
```
Help for functions can be obtained by typing:
`?functionName`
