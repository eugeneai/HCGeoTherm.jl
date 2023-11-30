using HCGT
using Test


function run()
    q0 = 33:0.2:40         # [mW/m^2] surface heat flow
    # q0 = 20:10:100         # [mW/m^2] surface heat flow
    GP = defaultGTInit(q0, true)
    dataf = userLoadCSV("./data/PT Ybileynaya_Gtherm.csv")
    answer = userComputeGeotherm(GP, dataf)
    userPlot(answer, GP.gfxDir,
             "geotherm.svg",
             "geotherm-chisquare.svg",
             "geotherm-opt.svg")
end

@testset "HCGT.jl" begin
    # Write your tests here.
end
