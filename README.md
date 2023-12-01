# HCGeoTherm.jl package

Geotherm calculation routines according to Dr. Derrek Hasterok and Dr. David Chapman

[![Build Status](https://github.com/eugeneai/HCGeoTherm.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/eugeneai/HCGeoTherm.jl/actions/workflows/CI.yml?query=branch%3Amaster)

At present it is a routines used in an application at https://gtherm.ru . Fill free improve package structure

### A basic considerations of usage

Here some examples of the library usage (actually taken from tests).  Data used in a test is published in a scientific journal and is subject of [CC BY-NC 4.0 DEED](https://creativecommons.org/licenses/by-nc/4.0/deed.en) license (c) [Anna Dymshits](https://www.researchgate.net/profile/Anna-Dymshits)

First, install the package (it pulls a lot of other):

```julia
Pkg.add("HCGeoTherm")
```

Then run example

```julia
using HCGeoTherm
using CSV
using DataFrames

# The CSV is expected to have columns T_C or T_K,   D_km or D_m or P_GPa, or P_kbar
# "P_GPa";"D_km";"T_C"
# 5,01351657010213;158,710903731105;1053,09689717222
# 6,27799615553237;197,151083128184;1189,07435289063
# 6,2812869402397;197,251122983287;1197,13815898029
# . . . . . . . . . . . . . . . . . . . . . . . . .
#

function loadCSV(fileName :: String, canonify::Bool = false) :: DataFrame
    pt = CSV.read(fileName, DataFrame, delim=';', decimal=',')
    if canonify
        pt |> canonifyDF    # Try to compute required columns
    end
    pt
end



q0 = 33:0.2:40                      # [mW/m^2] surface heat flow
GP = defaultGTInit(q0, true)        # true - switch on optimization
                                    # after the "true" You can add a folder, where put
                                    # output graphics, by default it is "/var/tmp"
df = loadCSV("data/PTdata.csv")     # Place here reference to CSV data file, see format below
                                    # The package should contain demo data,
                                    # or take it from github repo at
                                    # https://github.com/eugeneai/HCGTherm.jl/tree/master/test/data
println("Initial data frame")
println(df)
df.T_K = df.T_C .+ 273
answer = computeGeotherm(GP, df)
# print or plot results.            # We removed plotting to make the package cleaner and faster to start.
```

After probable successful run three SVG file will appear in ```/var/tmp```

(Here will go paper DOI's of corresponding publications)
