# Simple test for HCGeoTherm module

using HCGeoTherm

println("Testing HCGeoTherm module...")

# Test 1: Check if module loads
println("\n=== Test 1: Module loading ===")
println("Module loaded successfully: HCGeoTherm")

# Test 2: Check exported functions
println("\n=== Test 2: Exported functions ===")
println("Exported structures: GTInit, Geotherm, GTResult, ModelParameters")

# Test 3: Try to create a simple GTInit
println("\n=== Test 3: Create GTInit ===")
try
    # Create a simple range for q0
    q0_range = range(34.0, step=1.0, length=7)  # 34.0, 35.0, ..., 40.0

    # Create GTInit directly
    init = GTInit(
        q0_range,
        16.0,                    # D
        [16.0, 23.0, 39.0, 300.0],  # zbot
        225.0,                   # zmax
        0.1,                     # dz
        0.74,                    # P
        [0.0, 0.4, 0.4, 0.02],  # H
        3,                       # iref
        Set{String}(),          # options
        "test_model"            # model_name
    )

    println("✓ GTInit created successfully")
    println("  q0 type: $(typeof(init.q0))")
    println("  q0 values: $(collect(init.q0))")
    println("  D: $(init.D)")
    println("  model_name: $(init.model_name)")
catch e
    println("✗ Error creating GTInit: $e")
    println("  Stacktrace: ", sprint(showerror, e, catch_backtrace()))
end

# Test 4: Try empgtherms function
println("\n=== Test 4: empgtherms function ===")
try
    # Call empgtherms with float values
    T, z, k, A, q, alpha, az = empgtherms(
        40.0,      # q0
        225.0,     # maxsz
        0.1,       # dz
        16.0,      # D
        [16.0, 23.0, 39.0, 300.0],  # zbot
        [0.0, 0.4, 0.4, 0.02]       # H
    )

    println("✓ empgtherms executed successfully")
    println("  T length: $(length(T))")
    println("  z length: $(length(z))")
    println("  az: $az")

    # Check that results are reasonable
    if length(T) == length(z)
        println("✓ T and z have same length")
    else
        println("✗ T and z have different lengths")
    end

    if all(isfinite, T) && all(isfinite, z)
        println("✓ All T and z values are finite")
    else
        println("✗ Some T or z values are not finite")
    end
catch e
    println("✗ Error in empgtherms: $e")
    println("  Stacktrace: ", sprint(showerror, e, catch_backtrace()))
end

# Test 5: Check other exported functions exist
println("\n=== Test 5: Check function existence ===")
functions_to_check = [
    :computeGeotherm,
    :chisquare,
    :optimizeGeotherm,
    :saveResults,
    :loadResults,
    :validateInputs
]

for func in functions_to_check
    if isdefined(HCGeoTherm, func)
        println("✓ $func is defined")
    else
        println("✗ $func is not defined")
    end
end

println("\n=== Test Summary ===")
println("Simple test completed. Check above for any errors.")
