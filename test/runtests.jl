using HCGTherm
using Test
using DataFrames

function run()::GTResult
    q0 = 33:0.2:40         # [mW/m^2] surface heat flow
    GP = defaultGTInit(q0, true)
    df = userLoadCSV("data/PTdata.csv")
    println("Initial data frame")
    println(df)
    df.T_K = df.T_C .+ 273
    answer = userComputeGeotherm(GP, df)
    userPlot(answer, GP.gfxDir,
             "geotherm.svg",
             "geotherm-chisquare.svg",
             "geotherm-opt.svg")
    answer
end

@testset "compute without optimization" begin
    GP = defaultGTInit([35])
    P_GPa=[6.54]
    T_C=[1368]
    T_K=T_C .+ 273
    D_km=P_GPa |> depth
    df = DataFrame(P_GPa=P_GPa, T_C=T_C, T_K=T_K, D_km=D_km)
    answer = userComputeGeotherm(GP, df)
    T = answer.T
    lT = T[length(T)]
    @test lT>1000
    @test lT<2000
end


@testset "Complex all-in-one test" begin
    rc = run()
    @test true # Reached here
    dir = rc.ini.gfxDir * "/"
    println("Output dir:", dir)
    @test isfile(dir * "geotherm.svg")
    @test isfile(dir * "geotherm-chisquare.svg")
    @test isfile(dir * "geotherm-opt.svg")
end
