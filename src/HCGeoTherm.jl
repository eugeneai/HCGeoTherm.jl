module HCGeoTherm

using DataFrames
using Optim
using Formatting
using Interpolations

export
    GTInit, Geotherm, GTResult,
    defaultGTInit, depth, pressure,
    computeGeotherm,
    empgtherms, tccomp, thermcond, kcoef,
    empexpansivity, acoef,
    canonifyDF

include("EmpgTherm.jl")

struct GTInit  # Structure for model task representation
    q0   :: Any
        # StepRange{Int64, Int64}   # [mW/m^2] surface heat flow
    D    :: Float64                 # [km] thickness of upper crust
    zbot :: Vector{Int64}           # [km] base of lithospheric layers
    zmax :: Float64                 # [km] maximum depth of model
    dz   :: Float64                 # [km] depth step
    P    :: Float64                 # partition coefficient for upper crustal
                                    # heat production
    H    :: Vector{Float64}         # [uW/m^3] heat production of lithospheric
    iref :: Int64                   # index to reference heat flow
                                    # for elevation computation
    opt  :: Bool                    # Wether to apply optimization for q0
end

struct Geotherm
    T :: Vector{Float64}
    z :: Vector{Float64}
    label :: String
    q0:: Float64
    az:: Float64
end

struct GTResult
    ini :: GTInit
    GT :: Vector{Geotherm}
    GT_opt :: Union{Geotherm,Nothing}
    D :: Union{DataFrame,Nothing}
    max :: Union{DataFrame,Nothing}
end

function defaultGTInit(q0 = 34:1:40,
                       opt::Bool = false) :: GTInit
    GTInit(q0, 16, [16,23,39,300], 225,
           0.1, 0.74, [0,0.4,0.4,0.02], 3, opt)
end

function depth(p) # GPa -> km
    p .* 30.4 .+ 6.3
end

function pressure(d) # km -> GPa
    (d .- 6.3) ./ 30.4
end

function canonifyDF(pt::DataFrame)::DataFrame
    pt_n = names(pt)
    if "T_C" in pt_n
        TC = pt.T_C
        TK = TC .+ 273
    elseif "T_K" in pt_n
        TK = pt.T_K
        TC = TK .- 273
    end
    if "D_km" in pt_n
        Dkm = pt.D_km
        PGPa = Dkm |> pressure
    elseif "D_m" in pt_n
        Dkm = pt.D_m ./ 1000
        PGPa = Dkm |> pressure
    elseif "P_GPa" in pt_n
        PGPa = pt.P_GPa
        Dkm = PGPa |> depth
    elseif "P_kbar" in pt_n
        PGPa = pt.P_kbar ./ 10.0
        Dkm = PGPa |> depth
    end
    DataFrame(D_km=Dkm, P_GPa=PGPa, T_C=TC, T_K=TK)
end

function computeGeotherm(initParameters :: GTInit,
                         df :: Union{DataFrame,Nothing}=nothing) :: GTResult
    ini = initParameters
    if ! isnothing(df)
        maximumf = combine(df, [:D_km, :T_C, :T_K] .=> maximum)
    else
        maximumf = nothing
    end

    T = undef
    alpha_ = undef
    de=zeros(0)
    GTs = Vector{Geotherm}()
    for i = 1:length(ini.q0)
        # compute surface heat production
        A0 = (1 - ini.P) * ini.q0[i] / ini.D
        H = copy(ini.H)
        H[1] = A0

        # compute geotherm from emperical rather than long form
        _T,z,k,A,q,_alpha_,az = empgtherms(ini.q0[i],
                                        ini.zmax,
                                        ini.dz,
                                        ini.D,
                                        ini.zbot,
                                        H)
        if az>0.0
            label = format("{} ({})", ini.q0[i], az)
        else
            label = format("{}", ini.q0[i])
        end
        push!(GTs, Geotherm(_T, z, label, ini.q0[i], az))

        if T == undef
            T = _T
        else
            T = hcat(T,_T)
        end
        if alpha_ == undef
            alpha_ = _alpha_
        else
            alpha_ = hcat(alpha_,_alpha_)
        end

        # thermal elevation from paper (de) from emperical geotherms (dem)
        if length(ini.q0) > 1
            de = vcat(de, sum(T[:,i].*alpha_[:,i] - T[:,1].*alpha_[:,1])*ini.dz)
        end

    end
    answer = GTResult(ini, GTs, nothing, df, maximumf)

    if ini.opt
        (xs, ifu) = chisquare(answer)
        nxsb = minimum(xs)
        nxse = maximum(xs)
        nxs = nxsb:((nxse-nxsb)/100):nxse
        sp = convert(Float64, nxsb + (nxse-nxsb)/2)

        res = optimize(ifu,
                       convert(Float64, nxsb),
                       convert(Float64, nxse),
                       GoldenSection())

        miny = Optim.minimum(res)
        minx = Optim.minimizer(res)

        q0 = convert(Float64, minx)         # [mW/m^2] surface heat flow

        ai = answer.ini

        gpOpt = GTInit([q0], ai.D, ai.zbot, ai.zmax, ai.dz, ai.P, ai.H, ai.iref, false)

        answero = computeGeotherm(gpOpt, answer.D)

        GTResult(ini, GTs, answero.GT[1], df, maximumf)
    else
        answer
    end

end

function myInterpolate(xs, ys)
    b = minimum(xs)
    e = maximum(xs)
    dx = (e-b)/(length(xs)-1)
    ifu = interpolate(ys, BSpline(Cubic(Throw(OnGrid()))))
    function f(x)
        bb = x.-b
        v = bb./dx
        ifu(v.+1)
    end
    f
end

function chisquareGT(GT::Geotherm, D::DataFrame) :: Float64
    z = GT.z
    T = GT.T
    gti = myInterpolate(z,T)
    s = 0.0 :: Float64
    for row in eachrow(D)
        cT = (gti(row.D_km)-row.T_C)^2
        s = s + cT
    end
    s
end

function chisquare(result::GTResult) :: Any
    qs = result.ini.q0
    chis = Vector{Float64}()
    for (i, q) in enumerate(qs)
        chi = chisquareGT(result.GT[i], result.D)
        push!(chis, chi)
    end
    (qs, myInterpolate(qs, chis))
end

end # module
