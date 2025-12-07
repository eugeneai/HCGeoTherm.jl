# Test file to verify q0 as Float64 works correctly in HCGeoTherm

using HCGeoTherm

println("Testing q0 as Float64 in HCGeoTherm...")

function run_tests()
    passed = 0
    total = 0

    function check(condition, description)
        total += 1
        if condition
            println("✓ $description")
            passed += 1
        else
            println("✗ $description")
        end
    end

    println("\n=== Test 1: Default initialization with float range ===")
    try
        init = defaultGTInit()
        check(eltype(init.q0) == Float64, "q0 should be Float64")
        check(minimum(init.q0) == 34.0, "minimum q0 should be 34.0")
        check(maximum(init.q0) == 40.0, "maximum q0 should be 40.0")
        check(step(init.q0) == 1.0, "step should be 1.0")
    catch e
        println("Error in test 1: $e")
    end

    println("\n=== Test 2: Create GTInit with single float value ===")
    try
        init = createGTInit(40.5, 16.0, [16.0, 23.0, 39.0, 300.0], 225.0, 0.1, 0.74, [0.0, 0.4, 0.4, 0.02])
        check(eltype(init.q0) == Float64, "q0 should be Float64")
        check(minimum(init.q0) == 40.5, "minimum q0 should be 40.5")
        check(maximum(init.q0) == 40.5, "maximum q0 should be 40.5")
        check(length(collect(init.q0)) == 1, "should have 1 q0 value")
    catch e
        println("Error in test 2: $e")
    end

    println("\n=== Test 3: Create GTInit with float vector ===")
    try
        q0_values = [34.0, 35.5, 37.0, 38.5, 40.0]
        init = createGTInit(q0_values, 16.0, [16.0, 23.0, 39.0, 300.0], 225.0, 0.1, 0.74, [0.0, 0.4, 0.4, 0.02])
        check(eltype(init.q0) == Float64, "q0 should be Float64")
        check(minimum(init.q0) == 34.0, "minimum q0 should be 34.0")
        check(maximum(init.q0) == 40.0, "maximum q0 should be 40.0")
        check(abs(step(init.q0) - 1.5) < 1e-10, "step should be approximately 1.5")
    catch e
        println("Error in test 3: $e")
    end

    println("\n=== Test 4: Create GTInit with float StepRange ===")
    try
        init = createGTInit(34.0:0.5:40.0, 16.0, [16.0, 23.0, 39.0, 300.0], 225.0, 0.1, 0.74, [0.0, 0.4, 0.4, 0.02])
        check(eltype(init.q0) == Float64, "q0 should be Float64")
        check(minimum(init.q0) == 34.0, "minimum q0 should be 34.0")
        check(maximum(init.q0) == 40.0, "maximum q0 should be 40.0")
        check(step(init.q0) == 0.5, "step should be 0.5")
    catch e
        println("Error in test 4: $e")
    end

    println("\n=== Test 5: Geotherm calculation with float q0 ===")
    try
        init = createGTInit(34.0:0.5:36.0, 16.0, [16.0, 23.0, 39.0, 300.0], 225.0, 0.1, 0.74, [0.0, 0.4, 0.4, 0.02])
        result = computeGeotherm(init)

        check(length(result.GT) == 5, "should have 5 geotherms")
        check(all(gt -> gt.q0 isa Float64, result.GT), "all q0 values should be Float64")
        check(abs(result.GT[1].q0 - 34.0) < 1e-10, "first q0 should be 34.0")
        check(abs(result.GT[2].q0 - 34.5) < 1e-10, "second q0 should be 34.5")
        check(abs(result.GT[3].q0 - 35.0) < 1e-10, "third q0 should be 35.0")
    catch e
        println("Error in test 5: $e")
    end

    println("\n=== Test 6: Non-integer step size ===")
    try
        init = createGTInit(30.0:0.3:31.2, 16.0, [16.0, 23.0, 39.0, 300.0], 225.0, 0.1, 0.74, [0.0, 0.4, 0.4, 0.02])
        check(eltype(init.q0) == Float64, "q0 should be Float64")
        check(abs(minimum(init.q0) - 30.0) < 1e-10, "minimum q0 should be 30.0")
        check(abs(maximum(init.q0) - 31.2) < 1e-10, "maximum q0 should be 31.2")
        check(abs(step(init.q0) - 0.3) < 1e-10, "step should be 0.3")

        q0_collected = collect(init.q0)
        check(length(q0_collected) == 5, "should have 5 q0 values")
        check(abs(q0_collected[1] - 30.0) < 1e-10, "first value should be 30.0")
        check(abs(q0_collected[2] - 30.3) < 1e-10, "second value should be 30.3")
    catch e
        println("Error in test 6: $e")
    end

    println("\n=== Test 7: Direct GTInit constructor with float StepRange ===")
    try
        init = GTInit(
            34.0:0.25:35.0,
            16.0,
            [16.0, 23.0, 39.0, 300.0],
            225.0,
            0.1,
            0.74,
            [0.0, 0.4, 0.4, 0.02],
            3,
            Set{String}(),
            "test_model"
        )

        check(eltype(init.q0) == Float64, "q0 should be Float64")
        check(minimum(init.q0) == 34.0, "minimum q0 should be 34.0")
        check(maximum(init.q0) == 35.0, "maximum q0 should be 35.0")
        check(step(init.q0) == 0.25, "step should be 0.25")

        q0_collected = collect(init.q0)
        check(length(q0_collected) == 5, "should have 5 q0 values")
    catch e
        println("Error in test 7: $e")
    end

    println("\n=== Test 8: empgtherms with float q0 ===")
    try
        q0_val = 40.5
        T, z, k, A, q, alpha, az = empgtherms(
            q0_val,
            225.0,
            0.1,
            16.0,
            [16.0, 23.0, 39.0, 300.0],
            [0.0, 0.4, 0.4, 0.02]
        )

        check(q0_val isa Float64, "q0 should be Float64")
        check(length(T) == length(z), "T and z should have same length")
        check(all(isfinite, T), "all T values should be finite")
        check(all(isfinite, z), "all z values should be finite")
    catch e
        println("Error in test 8: $e")
    end

    println("\n=== Test 9: Edge cases for q0 ===")
    try
        # Very small step
        init1 = createGTInit(40.0:0.001:40.01, 16.0, [16.0, 23.0, 39.0, 300.0], 225.0, 0.1, 0.74, [0.0, 0.4, 0.4, 0.02])
        check(abs(step(init1.q0) - 0.001) < 1e-10, "step should be 0.001")

        # Single value with custom step
        init2 = createGTInit(42.5:0.5:42.5, 16.0, [16.0, 23.0, 39.0, 300.0], 225.0, 0.1, 0.74, [0.0, 0.4, 0.4, 0.02])
        q0_collected2 = collect(init2.q0)
        check(length(q0_collected2) == 1, "should have 1 q0 value")
        check(abs(q0_collected2[1] - 42.5) < 1e-10, "q0 value should be 42.5")
    catch e
        println("Error in test 9: $e")
    end

    println("\n" * "="^50)
    println("Summary: $passed/$total tests passed")

    return passed == total
end



# Run the tests
if abspath(PROGRAM_FILE) == @__FILE__
    success = run_tests()
    exit(success ? 0 : 1)
end
