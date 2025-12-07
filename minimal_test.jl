# Minimal test for empgtherms function

using HCGeoTherm

println("Testing empgtherms function...")

try
    # Simple test with minimal parameters
    T, z, k, A, q, alpha, az = empgtherms(
        40.0,      # q0
        50.0,      # maxsz - smaller for faster test
        1.0,       # dz - larger step for faster test
        16.0,      # D
        [16.0, 23.0, 39.0, 50.0],  # zbot - adjusted for maxsz
        [0.0, 0.4, 0.4, 0.02]       # H
    )

    println("SUCCESS: empgtherms executed")
    println("  T length: $(length(T))")
    println("  z length: $(length(z))")
    println("  First T: $(T[1])")
    println("  Last T: $(T[end])")
    println("  az: $az")

    # Basic validation
    if length(T) == length(z)
        println("✓ T and z have same length")
    else
        println("✗ T and z have different lengths")
    end

    if all(isfinite, T)
        println("✓ All T values are finite")
    else
        println("✗ Some T values are not finite")
    end

catch e
    println("ERROR: $e")
    println("Stacktrace:")
    showerror(stdout, e, catch_backtrace())
end

println("\nTest completed.")
