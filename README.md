# HCGTherm.jl package

Geotherm calculation routines according to Dr. Derrek Hasterok and Dr. David Chapman

[![Build Status](https://github.com/eugeneai/HCGTherm.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/eugeneai/HCGTherm.jl/actions/workflows/CI.yml?query=branch%3Amaster)

At present it is a routines used in an application at https://gtherm.ru . Fill free improve package structure

### A basic considerations of usage

Here some examples of the library usage (actually taken from tests).  Data used in a test is published in a scientific journal and is subject of [CC BY-NC 4.0 DEED](https://creativecommons.org/licenses/by-nc/4.0/deed.en) license (c) [Anna Dymshits](https://www.researchgate.net/profile/Anna-Dymshits)

First, install the package (it pulls a lot of other):

```julia
Pkg.add("HCGTherm")
```

Then run example

```julia
using HCGTherm

q0 = 33:0.2:40                      # [mW/m^2] surface heat flow
GP = defaultGTInit(q0, true)        # true - switch on optimization
                                    # after the "true" You can add a folder, where put
                                    # output graphics, by default it is "/var/tmp"
df = userLoadCSV("data/PTdata.csv") # Place here reference to CSV data file, see format below
                                    # The package should contain demo data,
                                    # or take it from github repo at
                                    # https://github.com/eugeneai/HCGTherm.jl/tree/master/test/data
println("Initial data frame")
println(df)
df.T_K = df.T_C .+ 273
answer = userComputeGeotherm(GP, df)
userPlot(answer, GP.gfxDir,         # Result graphics by default is at /var/tmp
         "geotherm.svg",
         "geotherm-chisquare.svg",
         "geotherm-opt.svg")
```

After probable successful run three SVG file will appear in ```/var/tmp```

(Here will go paper DOI's of corresponding publications)
