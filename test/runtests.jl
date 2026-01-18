using ScatteringKMatrix
using Test
using ScatteringKMatrix.StaticArrays


@testset "TwoBodyChannel" begin
    ch = TwoBodyChannel(1.1, 1.1)
    @test ch.m1 == 1.1 + 0.0im
    @test ch.m2 == 1.1 + 0.0im
    @test ch.L == 0

    # Test phase space calculation
    m = 3.0
    iρval = iρ(ch, m)
    @test iρval isa Complex
    @test abs(iρval) > 0

    # Test below threshold
    @test imag(iρ(ch, 1.0)) ≈ 0 atol = 1e-10
end

@testset "3x3 one-pole K-matrix" begin
    channels = SVector(
        TwoBodyChannel(1.1, 1.1),
        TwoBodyChannel(2.2, 2.2),
        TwoBodyChannel(1.3, 1.3),
    )
    MG = [(M=5.3, gs=[1.2, 0.48, 1.6])]
    K = KMatrix(MG)
    T = TMatrix(K, channels)

    @test nchannels(T) == 3
    @test npoles(T) == 1

    # Test amplitude calculation
    m = 5.0
    A = amplitude(T, m)
    @test size(A) == (3, 3)
    @test A ≈ transpose(A) # Check hermiticity
    # 
    @test isapprox(
        amplitude(T, 5.0),
        [
            0.198372+0.230421im 0.0793488+0.0921683im 0.264496+0.307228im
            0.0793488+0.0921683im 0.0317395+0.0368673im 0.105798+0.122891im
            0.264496+0.307228im 0.105798+0.122891im 0.352662+0.409637im
        ],
        atol=1e-6,
    )
end

@testset "2x2 two-pole K-matrix" begin
    channels = SVector(
        TwoBodyChannel(1.1, 1.1),
        TwoBodyChannel(1.3, 1.3),
    )
    MG = [
        (M=4.3, gs=[2.1, 0.0]),
        (M=6.3, gs=[0.0, 2.5]),
    ]
    K = KMatrix(MG)
    T = TMatrix(K, channels)

    @test nchannels(T) == 2
    @test npoles(T) == 2

    # Test decoupled case
    m = 5.0
    A = amplitude(T, m)
    @test abs(A[1, 2]) ≈ 0 atol = 1e-10
    @test abs(A[2, 1]) ≈ 0 atol = 1e-10
end

@testset "Production amplitude" begin
    channels = SVector(TwoBodyChannel(1.1, 1.1))
    MG = [
        (M=4.3, gs=[2.1]),
        (M=6.3, gs=[2.5]),
    ]
    K = KMatrix(MG)
    T = TMatrix(K, channels)

    # Test with default production couplings
    A = ProductionAmplitude(T, [1.2, 2.2])
    @test length(A.α_poles) == 2
    @test all(A.α_poles .== [1.2, 2.2])

    # Test with custom production couplings
    α = SVector(1.0, 2.0 * cis(π / 4))
    A_custom = ProductionAmplitude(T, α, SVector(0.0))

    m = 5.0
    amp = amplitude(A_custom, m)
    @test length(amp) == 1

    # Test individual pole contributions
    amp1 = production_pole(A_custom, m, 1)
    amp2 = production_pole(A_custom, m, 2)
    @test amp1 ≈ [-0.3068709069551295 + 0.06943242478375794im]
    @test amp2 ≈ [0.2807585353071957 + 0.17715198048327294im]
end

@testset "ChewMandestam function" begin
    m1, m2 = 1.1, 1.1
    m1sq, m2sq = m1^2, m2^2
    sth = (m1 + m2)^2
    io = 1e-6im

    # Test 1: Compare CM(s+io) and CM(s-io) below and above threshold
    # Below threshold
    s_below = 0.5 * sth  # Well below threshold
    cm_plus_below = ChewMandestam(s_below + io, m1sq, m2sq)
    cm_minus_below = ChewMandestam(s_below - io, m1sq, m2sq)
    diff_below = cm_plus_below - cm_minus_below
    # Difference should be proportional to io (may be smaller due to analytic continuation)
    @test abs(diff_below) < abs(io) * 2  # Should be on the order of io or smaller

    # Above threshold
    s_above = 1.5 * sth  # Well above threshold
    cm_plus_above = ChewMandestam(s_above + io, m1sq, m2sq)
    cm_minus_above = ChewMandestam(s_above - io, m1sq, m2sq)
    # Imaginary part should flip sign between the two
    @test imag(cm_plus_above) ≈ -imag(cm_minus_above) atol = 1e-10
    @test real(cm_plus_above) ≈ real(cm_minus_above) atol = 1e-10

    # Test 2: Compare imag(CM(s+io)) with imag(iρ(TwoBodyChannel, s+i0)) at multiple points
    ch_tbc = TwoBodyChannel(m1, m2)
    test_points = [0.5 * sth, 0.8 * sth, 1.2 * sth, 1.5 * sth, 2.0 * sth]

    for s_test in test_points
        m_test = sqrt(s_test)
        cm_val = ChewMandestam(s_test + io, m1sq, m2sq)
        irho_val = iρ(ch_tbc, m_test + 1e-6im)

        # Compare imaginary parts
        @test imag(cm_val) ≈ imag(irho_val) atol = 1e-5
    end

    # Test 3: Check normalization at high s - imaginary part should approach 1
    s_very_high = 10000.0
    cm_very_high = ChewMandestam(s_very_high + io, m1sq, m2sq)
    irho_very_high = iρ(ch_tbc, sqrt(s_very_high) + 1e-6im)
    @test abs(imag(cm_very_high)) ≈ 1.0 atol = 1e-3
    @test abs(imag(irho_very_high)) ≈ 1.0 atol = 1e-3
    @test imag(cm_very_high) ≈ imag(irho_very_high) atol = 1e-5
end

@testset "TwoBodyChewMandelstamChannel" begin
    ch = TwoBodyChewMandelstamChannel(1.1, 1.1)
    @test ch.m1 == 1.1 + 0.0im
    @test ch.m2 == 1.1 + 0.0im
    @test ch.L == 0
    @test threshold(ch) ≈ 2.2

    # Test phase space calculation
    m = 3.0
    iρval = iρ(ch, m)
    @test iρval isa Complex
    @test abs(iρval) > 0

    # Test below threshold
    @test abs(imag(iρ(ch, 1.0))) < 1e-10

    # Test error handling for L != 0
    ch_l1 = TwoBodyChewMandelstamChannel(1.1, 1.1; L=1)
    @test_throws ErrorException iρ(ch_l1, 3.0)
end
