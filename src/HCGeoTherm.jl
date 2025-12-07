# HCGeoTherm.jl - Geotherm Modeling Module
# Enhanced version with improved structure, documentation, and error handling

module HCGeoTherm

using DataFrames
using Optim
using Format
using Interpolations
import Logging
using Statistics
using LinearAlgebra

export
    # Structures
    GTInit, Geotherm, GTResult, ModelParameters,

    # Initialization functions
    defaultGTInit, createModelParameters,

    # Core computation functions
    computeGeotherm, computeGeothermSeries,

    # Utility functions
    depth, pressure, canonifyDF,

    # Analysis functions
    chisquare, chisquareGT, optimizeGeotherm,

    # Export functions
    saveResults, loadResults,

    # Validation functions
    validateInputs, validateResults

# ============================================================================
# DATA STRUCTURES
# ============================================================================

"""
    GTInit

Structure representing initialization parameters for geotherm calculation.

# Fields
- `q0::StepRange{Int64, Int64}`: Surface heat flow range [mW/m²]
- `D::Float64`: Thickness of upper crust [km]
- `zbot::Vector{Float64}`: Depth to layer bases [km]
- `zmax::Float64`: Maximum depth of model [km]
- `dz::Float64`: Depth step [km]
- `P::Float64`: Partition coefficient for upper crustal heat production
- `H::Vector{Float64}`: Heat production of lithospheric layers [μW/m³]
- `iref::Int64`: Index to reference heat flow for elevation computation
- `options::Set{String}`: Optional calculation variants (e.g., {"opt", "misfits"})
- `model_name::String`: Name of the model configuration
"""
struct GTInit
    q0::StepRange{Int64, Int64}
    D::Float64
    zbot::Vector{Float64}
    zmax::Float64
    dz::Float64
    P::Float64
    H::Vector{Float64}
    iref::Int64
    options::Set{String}
    model_name::String

    function GTInit(q0, D, zbot, zmax, dz, P, H, iref, options, model_name="default")
        # Validate inputs
        @assert length(zbot) == length(H) "zbot and H must have same length"
        @assert all(zbot .> 0) "zbot must be positive"
        @assert all(H .>= 0) "H must be non-negative"
        @assert D > 0 "D must be positive"
        @assert zmax > 0 "zmax must be positive"
        @assert dz > 0 "dz must be positive"
        @assert 0 <= P <= 1 "P must be between 0 and 1"
        @assert iref >= 1 "iref must be at least 1"

        new(q0, D, zbot, zmax, dz, P, H, iref, options, model_name)
    end
end

"""
    Geotherm

Structure representing a single geotherm calculation result.

# Fields
- `T::Vector{Float64}`: Temperature profile [°C]
- `z::Vector{Float64}`: Depth profile [km]
- `label::String`: Label for the geotherm
- `q0::Float64`: Surface heat flow [mW/m²]
- `az::Float64`: Adiabat depth [km] (negative if not reached)
- `metadata::Dict{String, Any}`: Additional metadata
"""
struct Geotherm
    T::Vector{Float64}
    z::Vector{Float64}
    label::String
    q0::Float64
    az::Float64
    metadata::Dict{String, Any}

    function Geotherm(T, z, label, q0, az, metadata=Dict{String, Any}())
        @assert length(T) == length(z) "T and z must have same length"
        @assert all(z .>= 0) "z must be non-negative"
        new(T, z, label, q0, az, metadata)
    end
end

"""
    GTResult

Structure representing complete geotherm calculation results.

# Fields
- `ini::GTInit`: Initialization parameters
- `GT::Vector{Geotherm}`: Calculated geotherms
- `GT_opt::Union{Geotherm, Nothing}`: Optimized geotherm (if calculated)
- `D::Union{DataFrame, Nothing}`: Input data
- `max::Union{DataFrame, Nothing}`: Maximum values from input data
- `statistics::Dict{String, Any}`: Calculation statistics
- `timestamp::DateTime`: Calculation timestamp
"""
struct GTResult
    ini::GTInit
    GT::Vector{Geotherm}
    GT_opt::Union{Geotherm, Nothing}
    D::Union{DataFrame, Nothing}
    max::Union{DataFrame, Nothing}
    statistics::Dict{String, Any}
    timestamp::DateTime

    function GTResult(ini, GT, GT_opt, D, max, statistics=Dict{String, Any}(), timestamp=now())
        new(ini, GT, GT_opt, D, max, statistics, timestamp)
    end
end

"""
    ModelParameters

Structure for advanced model parameters.

# Fields
- `thermal_conductivity_model::String`: Thermal conductivity model name
- `adiabat_gradient::Float64`: Adiabatic gradient [K/km]
- `surface_temperature::Float64`: Surface temperature [°C]
- `adiabat_surface_temp::Float64`: Adiabatic temperature at surface [°C]
- `crust_density::Float64`: Crust density [kg/m³]
- `mantle_density::Float64`: Mantle density [kg/m³]
- `gravity::Float64`: Gravity acceleration [m/s²]
"""
struct ModelParameters
    thermal_conductivity_model::String
    adiabat_gradient::Float64
    surface_temperature::Float64
    adiabat_surface_temp::Float64
    crust_density::Float64
    mantle_density::Float64
    gravity::Float64

    function ModelParameters(;
        thermal_conductivity_model="hasterok_chapman_2011",
        adiabat_gradient=0.3,
        surface_temperature=0.0,
        adiabat_surface_temp=1300.0,
        crust_density=2850.0,
        mantle_density=3340.0,
        gravity=10.0
    )
        @assert adiabat_gradient > 0 "adiabat_gradient must be positive"
        @assert crust_density > 0 "crust_density must be positive"
        @assert mantle_density > 0 "mantle_density must be positive"
        @assert gravity > 0 "gravity must be positive"

        new(
            thermal_conductivity_model,
            adiabat_gradient,
            surface_temperature,
            adiabat_surface_temp,
            crust_density,
            mantle_density,
            gravity
        )
    end
end

# ============================================================================
# INITIALIZATION FUNCTIONS
# ============================================================================

"""
    defaultGTInit(q0=34:1:40, options::Set{String}=Set()) :: GTInit

Create default GTInit structure with standard parameters.

# Arguments
- `q0`: Surface heat flow range [mW/m²]
- `options`: Calculation options

# Returns
- `GTInit`: Default initialization parameters
"""
function defaultGTInit(q0::StepRange{Int64, Int64}=34:1:40,
                       options::Set{String}=Set()) :: GTInit
    GTInit(
        q0,
        16.0,                    # D: Thickness of upper crust [km]
        [16.0, 23.0, 39.0, 300.0],  # zbot: Layer boundaries [km]
        225.0,                   # zmax: Maximum depth [km]
        0.1,                     # dz: Depth step [km]
        0.74,                    # P: Partition coefficient
        [0.0, 0.4, 0.4, 0.02],  # H: Heat production [μW/m³]
        3,                       # iref: Reference index
        options,
        "default_model"
    )
end

"""
    createModelParameters(; kwargs...) :: ModelParameters

Create model parameters with optional custom values.

# Keyword Arguments
- `thermal_conductivity_model`: Thermal conductivity model name
- `adiabat_gradient`: Adiabatic gradient [K/km]
- `surface_temperature`: Surface temperature [°C]
- `adiabat_surface_temp`: Adiabatic temperature at surface [°C]
- `crust_density`: Crust density [kg/m³]
- `mantle_density`: Mantle density [kg/m³]
- `gravity`: Gravity acceleration [m/s²]

# Returns
- `ModelParameters`: Model parameters structure
"""
createModelParameters(; kwargs...) = ModelParameters(; kwargs...)

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

"""
    depth(p) :: Vector{Float64}

Convert pressure to depth using Birch's law approximation.

# Arguments
- `p::Vector{Float64}`: Pressure [GPa]

# Returns
- `Vector{Float64}`: Depth [km]
"""
function depth(p::Vector{Float64}) :: Vector{Float64}
    @assert all(p .>= 0) "Pressure must be non-negative"
    p .* 30.4 .+ 6.3
end

depth(p::Float64) :: Float64 = p * 30.4 + 6.3

"""
    pressure(d) :: Vector{Float64}

Convert depth to pressure using Birch's law approximation.

# Arguments
- `d::Vector{Float64}`: Depth [km]

# Returns
- `Vector{Float64}`: Pressure [GPa]
"""
function pressure(d::Vector{Float64}) :: Vector{Float64}
    @assert all(d .>= 0) "Depth must be non-negative"
    (d .- 6.3) ./ 30.4
end

pressure(d::Float64) :: Float64 = (d - 6.3) / 30.4

"""
    canonifyDF(pt::DataFrame) :: DataFrame

Convert various pressure-temperature data formats to canonical format.

# Arguments
- `pt::DataFrame`: Input pressure-temperature data

# Returns
- `DataFrame`: Canonical format with columns: D_km, P_GPa, T_C, T_K

# Supported input formats:
- Temperature: T_C (°C) or T_K (K)
- Depth/Pressure: D_km (km), D_m (m), P_GPa (GPa), P_kbar (kbar)
"""
function canonifyDF(pt::DataFrame) :: DataFrame
    pt_n = names(pt)

    # Validate required columns
    has_temp = "T_C" in pt_n || "T_K" in pt_n
    has_depth_pressure = "D_km" in pt_n || "D_m" in pt_n || "P_GPa" in pt_n || "P_kbar" in pt_n

    @assert has_temp "DataFrame must contain either T_C or T_K column"
    @assert has_depth_pressure "DataFrame must contain depth or pressure column"

    # Convert temperature
    if "T_C" in pt_n
        TC = pt.T_C
        TK = TC .+ 273.15
    else  # "T_K" in pt_n
        TK = pt.T_K
        TC = TK .- 273.15
    end

    # Convert depth/pressure
    if "D_km" in pt_n
        Dkm = pt.D_km
        PGPa = pressure.(Dkm)
    elseif "D_m" in pt_n
        Dkm = pt.D_m ./ 1000.0
        PGPa = pressure.(Dkm)
    elseif "P_GPa" in pt_n
        PGPa = pt.P_GPa
        Dkm = depth.(PGPa)
    else  # "P_kbar" in pt_n
        PGPa = pt.P_kbar ./ 10.0
        Dkm = depth.(PGPa)
    end

    # Create canonical DataFrame
    DataFrame(D_km=Dkm, P_GPa=PGPa, T_C=TC, T_K=TK)
end

# ============================================================================
# CORE COMPUTATION FUNCTIONS
# ============================================================================

"""
    computeGeotherm(initParameters::GTInit, df::Union{DataFrame, Nothing}=nothing;
                    model_params::ModelParameters=ModelParameters()) :: Dict{String, GTResult}

Compute geotherms for given initialization parameters.

# Arguments
- `initParameters::GTInit`: Initialization parameters
- `df::Union{DataFrame, Nothing}`: Input data for optimization (optional)
- `model_params::ModelParameters`: Model physical parameters

# Returns
- `Dict{String, GTResult}`: Dictionary of results with keys: "series", "optimize", "misfits"

# Throws
- `ArgumentError` if inputs are invalid
"""
function computeGeotherm(initParameters::GTInit,
                         df::Union{DataFrame, Nothing}=nothing;
                         model_params::ModelParameters=ModelParameters()) :: Dict{String, GTResult}

    # Validate inputs
    validateInputs(initParameters, df)

    ini = initParameters
    statistics = Dict{String, Any}()

    # Prepare maximum values from input data
    maximumf = isnothing(df) ? nothing : combine(df, [:D_km, :T_C, :T_K] .=> maximum)

    # Calculate geotherm series
    GTs, calculation_stats = computeGeothermSeries(ini, model_params)
    merge!(statistics, calculation_stats)

    # Create base result
    answer = GTResult(ini, GTs, nothing, df, maximumf, statistics)

    answers = Dict{String, GTResult}("series" => answer)

    # Optimization if requested
    if "optimize" in ini.options || "misfits" in ini.options
        try
            optimized_result = optimizeGeotherm(answer, model_params)
            merge!(answers, optimized_result)
        catch e
            @warn "Optimization failed: $(e)"
        end
    end

    return answers
end

"""
    computeGeothermSeries(ini::GTInit, model_params::ModelParameters) :: Tuple{Vector{Geotherm}, Dict{String, Any}}

Compute series of geotherms for different surface heat flows.

# Arguments
- `ini::GTInit`: Initialization parameters
- `model_params::ModelParameters`: Model physical parameters

# Returns
- `Tuple{Vector{Geotherm}, Dict{String, Any}}`: Geotherms and calculation statistics
"""
function computeGeothermSeries(ini::GTInit, model_params::ModelParameters) :: Tuple{Vector{Geotherm}, Dict{String, Any}}
    GTs = Vector{Geotherm}()
    statistics = Dict{String, Any}()

    start_time = time()
    adiabat_count = 0

    for (i, q0_val) in enumerate(ini.q0)
        # Compute surface heat production
        A0 = (1 - ini.P) * q0_val / ini.D
        H = copy(ini.H)
        H[1] = A0

        # Compute geotherm
        T, z, _, _, _, _, az = empgtherms(
            Float64(q0_val), ini.zmax, ini.dz, ini.D,
            ini.zbot, H, model_params
        )

        # Create label
        if az > 0.0
            label = format("{:.2f} ({:.2f})", q0_val, az)
            adiabat_count += 1
        else
            label = format("{:.2f}", q0_val)
        end

        # Create geotherm with metadata
        metadata = Dict{String, Any}(
            "index" => i,
            "adiabat_reached" => az > 0.0,
            "calculation_time" => time() - start_time
        )

        push!(GTs, Geotherm(T, z, label, Float64(q0_val), az, metadata))
    end

    # Update statistics
    statistics["total_calculation_time"] = time() - start_time
    statistics["adiabat_count"] = adiabat_count
    statistics["total_geotherms"] = length(GTs)
    statistics["heat_flow_range"] = (minimum(ini.q0), maximum(ini.q0))

    return GTs, statistics
end

# ============================================================================
# OPTIMIZATION FUNCTIONS
# ============================================================================

"""
    chisquareGT(GT::Geotherm, D::DataFrame; sigmaT=20.0, sigmaP=0.3) :: Float64

Calculate chi-square misfit between geotherm and data.

# Arguments
- `GT::Geotherm`: Geotherm to evaluate
- `D::DataFrame`: Data for comparison
- `sigmaT`: Temperature uncertainty [K]
- `sigmaP`: Pressure uncertainty [GPa]

# Returns
- `Float64`: Chi-square misfit value
"""
function chisquareGT(GT::Geotherm, D::DataFrame; sigmaT=20.0, sigmaP=0.3) :: Float64
    @assert !isempty(D) "DataFrame cannot be empty"
    @assert sigmaT > 0 "sigmaT must be positive"
    @assert sigmaP > 0 "sigmaP must be positive"

    z = GT.z
    T = GT.T

    # Create interpolators
    T_interp = linear_interpolation(z, T; extrapolation_bc=Line())
    z_interp = linear_interpolation(T, z; extrapolation_bc=Line())

    dT = 0.0
    dP = 0.0

    for row in eachrow(D)
        # Temperature misfit
        T_pred = T_interp(row.D_km)
        dT += (T_pred - row.T_C)^2

        # Pressure misfit
        z_pred = z_interp(row.T_C)
        P_pred = pressure(z_pred)
        P_actual = pressure(row.D_km)
        dP += (P_pred - P_actual)^2
    end

    # Normalized chi-square
    n = nrow(D)
    chi2 = (dP / (sigmaP^2) + dT / (sigmaT^2)) / n

    return sqrt(chi2)
end

"""
    chisquare(result::GTResult; sigmaT=20.0, sigmaP=0.3) :: Tuple{Vector{Float64}, Function}

Calculate chi-square values for all geotherms in result.

# Arguments
- `result::GTResult`: Geotherm results
- `sigmaT`: Temperature uncertainty [K]
- `sigmaP`: Pressure uncertainty [GPa]

# Returns
- `Tuple{Vector{Float64}, Function}`: Heat flow values and interpolated chi-square function
"""
function chisquare(result::GTResult; sigmaT=20.0, sigmaP=0.3) :: Tuple{Vector{Float64}, Function}
    @assert !isnothing(result.D) "Result must contain data for chi-square calculation"

    qs = [gt.q0 for gt in result.GT]
    chis = Vector{Float64}()

    for gt in result.GT
        chi = chisquareGT(gt, result.D; sigmaT=sigmaT, sigmaP=sigmaP)
        push!(chis, chi)
    end

    # Create interpolating function
    chi_interp = linear_interpolation(qs, chis; extrapolation_bc=Line())

    return (qs, chi_interp)
end

"""
    optimizeGeotherm(result::GTResult, model_params::ModelParameters) :: Dict{String, GTResult}

Optimize geotherm to minimize chi-square misfit.

# Arguments
- `result::GTResult`: Initial geotherm results
- `model_params::ModelParameters`: Model physical parameters

# Returns
- `Dict{String, GTResult}`: Optimized results
"""
function optimizeGeotherm(result::GTResult, model_params::ModelParameters) :: Dict{String, GTResult}
    @assert !isnothing(result.D) "Data required for optimization"

    answers = Dict{String, GTResult}()

    # Get chi-square function
    qs, chi_func = chisquare(result)

    # Find minimum using Golden Section search
    q_min = minimum(qs)
    q_max = maximum(qs)

    res = optimize(chi_func, Float64(q_min), Float64(q_max), GoldenSection())

    if !Optim.converged(res)
        @warn "Optimization did not converge"
        return answers
    end

    optimal_q = Optim.minimizer(res)
    min_chi = Optim.minimum(res)

    # Create optimized initialization
    ai = result.ini
    gpOpt = GTInit(
        [round(Int, optimal_q)],
        ai.D, ai.zbot, ai.zmax, ai.dz,
        ai.P, ai.H, ai.iref,
        Set{String}(),
        "optimized_$(optimal_q)"
    )

    # Compute optimized geotherm
    opt_results = computeGeotherm(gpOpt, result.D; model_params=model_params)
    opt_result = opt_results["series"]

    # Update statistics
    stats = copy(result.statistics)
    stats["optimal_heat_flow"] = optimal_q
    stats["minimum_chi_square"] = min_chi
    stats["optimization_converged"] = Optim.converged(res)

    # Create result with optimized geotherm
    answer_opt = GTResult(
        result.ini,
        result.GT,
        opt_result.GT[1],
        result.D,
        result.max,
        stats,
        now()
    )

    answers["optimize"] = answer_opt

    # Calculate misfit bounds if requested
    if "misfits" in result.ini.options
        chi_std = sqrt(min_chi)
        q_low = optimal_q - chi_std
        q_high = optimal_q + chi_std

        gpMisfit = GTInit(
            [round(Int, q_low), round(Int, optimal_q), round(Int, q_high)],
            ai.D, ai.zbot, ai.zmax, ai.dz,
            ai.P, ai.H, ai.iref,
            Set{String}(),
            "misfit_bounds"
        )

        misfit_results = computeGeotherm(gpMisfit, result.D; model_params=model_params)
        misfit_result = misfit_results["series"]

        misfit_stats = copy(stats)
        misfit_stats["misfit_lower_bound"] = q_low
        misfit_stats["misfit_upper_bound"] = q_high

        answer_misfit = GTResult(
            result.ini,
            misfit_result.GT,
            opt_result.GT[1],
            result.D,
            result.max,
            misfit_stats,
            now()
        )

        answers["misfits"] = answer_misfit
    end

    return answers
end

# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

"""
    validateInputs(ini::GTInit, df::Union{DataFrame, Nothing})

Validate input parameters for geotherm calculation.

# Arguments
- `ini::GTInit`: Initialization parameters
- `df::Union{DataFrame, Nothing}`: Input data

# Throws
- `ArgumentError` if validation fails
"""
function validateInputs(ini::GTInit, df::Union{DataFrame, Nothing})
    # Validate GTInit
    @assert ini.zmax > maximum(ini.zbot) "zmax must be greater than all zbot values"
    @assert ini.dz < ini.zmax "dz must be smaller than zmax"
    @assert ini.iref <= length(ini.q0) "iref must be within q0 range"

    # Validate DataFrame if provided
    if !isnothing(df)
        @assert !isempty(df) "DataFrame cannot be empty"
        required_cols = ["D_km", "T_C"]
        for col in required_cols
            @assert col in names(df) "DataFrame must contain $col column"
        end

        @assert all(df.D_km .>= 0) "Depth must be non-negative"
        @assert all(df.T_C .>= -273.15) "Temperature must be above absolute zero"
    end
end

"""
    validateResults(result::GTResult) :: Bool

Validate geotherm calculation results.

# Arguments
- `result::GTResult`: Calculation results

# Returns
- `Bool`: True if results are valid
"""
function validateResults(result::GTResult) :: Bool
    try
        # Check geotherms
        @assert !isempty(result.GT) "No geotherms calculated"

        for gt in result.GT
            @assert length(gt.T) == length(gt.z) "Temperature and depth arrays must have same length"
            @assert all(gt.z .>= 0) "Depth must be non-negative"
            @assert all(gt.T .>= -273.15) "Temperature must be above absolute zero"
        end

        # Check optimized geotherm if present
        if !isnothing(result.GT_opt)
            @assert length(result.GT_opt.T) == length(result.GT_opt.z) "Optimized geotherm arrays must have same length"
        end

        return true
    catch e
        @warn "Validation failed: $(e)"
        return false
    end
end

# ============================================================================
# IMPORT/EXPORT FUNCTIONS
# ============================================================================

"""
    saveResults(result::GTResult, filepath::String)

Save geotherm results to JSON file.

# Arguments
- `result::GTResult`: Results to save
- `filepath::String`: Path to save file
"""
function saveResults(result::GTResult, filepath::String)
    # This would be implemented with JSON3 or similar package
    @warn "saveResults not implemented - requires JSON serialization"
end

"""
    loadResults(filepath::String) :: GTResult

Load geotherm results from JSON file.

# Arguments
- `filepath::String`: Path to load file

# Returns
- `GTResult`: Loaded results
"""
function loadResults(filepath::String) :: GTResult
    # This would be implemented with JSON3 or similar package
    error("loadResults not implemented - requires JSON deserialization")
end

# ============================================================================
# MODULE INITIALIZATION
# ============================================================================

# Include submodules
include("EmpgTherm.jl")

# Initialize module
function __init__()
    @info "HCGeoTherm module initialized"
    @info "Version: 2.0.0"
    @info "Enhanced with better structure, documentation, and error handling"
end

end # module HCGeoTherm
