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

    # Test 4: Check real part normalization - should be zero at threshold
    io_small = 1e-9im
    cm_at_threshold = ChewMandestam(sth + io_small, m1sq, m2sq)
    @test abs(cm_at_threshold) < 1.5e-5
end

@testset "TwoBodyChewMandelstamChannel" begin
    ch = TwoBodyChewMandelstamChannel(1.1, 1.1)
    @test TwoBodyChannel(ch) == TwoBodyChannel(1.1, 1.1)
    @test TwoBodyChannel(TwoBodyChewMandelstamChannel(0.14, 1.87)) ==
          TwoBodyChannel(0.14, 1.87)
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

@testset "K-matrix with TwoBodyChewMandelstamChannel" begin
    # Create channels using Chew-Mandelstam
    channels_cm = SVector(
        TwoBodyChewMandelstamChannel(1.1, 1.1),
        TwoBodyChewMandelstamChannel(2.2, 2.2),
    )

    # Create K-matrix and T-matrix
    MG = [(M=5.3, gs=[1.2, 0.48])]
    K = KMatrix(MG)
    T = TMatrix(K, channels_cm)

    @test nchannels(T) == 2
    @test npoles(T) == 1

    # Test amplitude calculation
    m = 5.0
    A = amplitude(T, m)
    @test size(A) == (2, 2)
    @test A ≈ transpose(A) # Check hermiticity
    #
    @test isapprox(
        amplitude(T, 5.0),
        [
            0.30102044275247963-0.09751950975706329im 0.12040817710099183-0.03900780390282532im
            0.12040817710099186-0.03900780390282533im 0.04816327084039675-0.015603121561130133im
        ],
        atol=1e-6,
    )
end

@testset "Square-root phase-space analytic continuation" begin
    for (m1, m2) in ((1.1, 1.1), (0.14, 1.87))
        ch = TwoBodyChannel(m1, m2)
        first_sheet = ContinuationChannel(ch)
        second_sheet = ContinuationChannel(ch, 2)
        glued_sheet = ContinuationChannel(ch, 12)
        cut_sheet = ContinuationChannel(ch, -90)
        th = threshold(ch)

        @test threshold(first_sheet) == th
        @test_throws ArgumentError ContinuationChannel(ch, 3)
        @test_throws ErrorException iρ(
            ContinuationChannel(TwoBodyChannel(m1, m2; L=1), 2),
            th + 1,
        )

        for m in (th + 0.3, th + 1.0, 2 * th)
            ε = 1e-7
            @test iρ(second_sheet, m) ≈ -iρ(first_sheet, m)
            @test iρ(first_sheet, m + im * ε) ≈ iρ(second_sheet, m - im * ε) atol = 1e-5
            @test abs(iρ(glued_sheet, m + im * ε) - iρ(glued_sheet, m - im * ε)) < 1e-5
            @test abs(iρ(cut_sheet, m + im * ε) - iρ(cut_sheet, m - im * ε)) < 1e-5
        end

        # The legacy square-root prescription is precisely the downward-cut
        # continuation, including away from the real axis.
        for m in (
            th + 0.5 + 0.2im,
            th + 0.5 - 0.2im,
            0.7 * th + 0.2im,
            0.7 * th - 0.2im,
        )
            @test iρ(cut_sheet, m) ≈ iρ(ch, m)
        end

        for α in (0, -π / 6, -π / 2, -π, π / 6)
            angled = AngledCutChannel(ch, α)
            mode = α == 0 ? 1 : α == -π / 2 ? -90 : α == -π ? 12 : nothing
            mode === nothing || @test iρ(angled, th + 0.5 - 0.1im) ==
                                      iρ(ContinuationChannel(ch, mode), th + 0.5 - 0.1im)
        end
    end

    channels_sqrt = SVector(
        TwoBodyChannel(1.5, 1.5),
        TwoBodyChannel(0.5, 0.5),
        TwoBodyChannel(1.0, 1.0),
    )
    continued = continue_channels(channels_sqrt, 2.5)
    @test getproperty.(continued, :mode) == SVector(1, 12, 12)
    @test eltype(continued) <: ContinuationChannel

    K = KMatrix([(M=3.5, gs=[1.0, 0.5, 0.2])])
    physical = TMatrix(K, map(ch -> ContinuationChannel(ch, 1), channels_sqrt))
    unphysical = TMatrix(K, continued)
    @test amplitude(unphysical, 2.4 + 0.1im) ≈ amplitude(physical, 2.4 + 0.1im)
    @test !(amplitude(unphysical, 2.4 - 0.1im) ≈ amplitude(physical, 2.4 - 0.1im))
    @test isfinite(detD(unphysical, 2.4 - 0.1im))
end

@testset "Chew-Mandelstam analytic continuation" begin
    # Sheet II differs from sheet I by the right-hand-cut discontinuity.
    # Modes 12 and -90 then select where that second sheet is used.
    for (m1, m2) in ((1.1, 1.1), (0.14, 1.87))
        ch = TwoBodyChewMandelstamChannel(m1, m2)
        first_sheet = ContinuationChannel(ch)
        second_sheet = ContinuationChannel(ch, 2)
        glued_sheet = ContinuationChannel(ch, 12)
        cut_sheet = ContinuationChannel(ch, -90)
        th = threshold(ch)

        @test first_sheet.mode == 1
        @test threshold(first_sheet) == th
        @test_throws ArgumentError ContinuationChannel(ch, 3)
        @test_throws ArgumentError ContinuationChannel(ch, 90)
        @test_throws ErrorException iρ(
            ContinuationChannel(TwoBodyChewMandelstamChannel(m1, m2; L=1), 2),
            th + 1,
        )

        for m in (th + 0.3, th + 1.0, 2 * th)
            ε = 1e-7
            jump = iρ(second_sheet, m) - iρ(first_sheet, m)
            disc = iρ(ch, m + im * ε) - iρ(ch, m - im * ε)
            @test jump ≈ disc atol = 1e-5
            @test real(jump) ≈ 0 atol = 1e-8
            @test imag(jump) > 0

            # Continuation through the cut: II just below matches I just above.
            @test iρ(first_sheet, m + im * ε) ≈ iρ(second_sheet, m - im * ε) atol = 1e-4
            @test abs(iρ(glued_sheet, m + im * ε) - iρ(glued_sheet, m - im * ε)) < 1e-4
            @test abs(iρ(cut_sheet, m + im * ε) - iρ(cut_sheet, m - im * ε)) < 1e-4

            # The real axis and the upper half-plane stay on sheet I,
            # except for the fixed second sheet.
            @test iρ(first_sheet, m) == iρ(ch, m)
            @test iρ(glued_sheet, m) == iρ(ch, m)
            @test iρ(cut_sheet, m) == iρ(ch, m)
            @test !(iρ(second_sheet, m) ≈ iρ(ch, m))
            @test iρ(glued_sheet, m + 0.2im) == iρ(first_sheet, m + 0.2im)
            @test iρ(cut_sheet, m + 0.2im) == iρ(first_sheet, m + 0.2im)
            @test iρ(glued_sheet, m - 0.2im) == iρ(second_sheet, m - 0.2im)
            @test iρ(cut_sheet, m - 0.2im) == iρ(second_sheet, m - 0.2im)
        end

        # Below threshold only mode 12 is on sheet II. Mode -90 keeps the
        # physical cut, so it still tracks sheet I.
        m_below = 0.7 * th
        @test iρ(glued_sheet, m_below + 0.2im) == iρ(first_sheet, m_below + 0.2im)
        @test iρ(glued_sheet, m_below - 0.2im) == iρ(second_sheet, m_below - 0.2im)
        @test iρ(cut_sheet, m_below - 0.2im) == iρ(first_sheet, m_below - 0.2im)
        @test abs(iρ(glued_sheet, m_below + im * 1e-4) - iρ(glued_sheet, m_below - im * 1e-4)) > 0.1
        @test abs(iρ(cut_sheet, m_below + im * 1e-4) - iρ(cut_sheet, m_below - im * 1e-4)) < 1e-3
        # The branch point itself stays on sheet I.
        @test iρ(cut_sheet, th - im * 1e-4) == iρ(first_sheet, th - im * 1e-4)
    end

    channels_cm = SVector(
        TwoBodyChewMandelstamChannel(1.5, 1.5),
        TwoBodyChewMandelstamChannel(0.5, 0.5),
        TwoBodyChewMandelstamChannel(1.0, 1.0),
    )
    # Thresholds are 3, 1, 2. Selection follows the threshold, not storage order.
    @test getproperty.(continue_channels(channels_cm, 0.2), :mode) == SVector(1, 1, 1)
    @test getproperty.(continue_channels(channels_cm, 1.0), :mode) == SVector(1, 12, 1)
    @test getproperty.(continue_channels(channels_cm, 2.5), :mode) == SVector(1, 12, 12)
    @test getproperty.(continue_channels(channels_cm, 3.0), :mode) == SVector(12, 12, 12)
    @test getproperty.(continue_channels(channels_cm, 2.5; mode=2), :mode) == SVector(1, 2, 2)
    @test getproperty.(continue_channels(channels_cm, 2.5; mode=-90), :mode) == SVector(1, -90, -90)
    @test_throws ArgumentError continue_channels(SVector(1.0), 3.0)
    @test_throws ArgumentError continue_channels(channels_cm, 2.5; mode=3)

    K = KMatrix([(M=3.5, gs=[1.0, 0.5, 0.2])])
    physical = TMatrix(K, channels_cm)
    continued = TMatrix(K, continue_channels(channels_cm, 2.5))
    # Reference mass 2.5 continues the channels with thresholds 1 and 2.
    # Above the real axis, and on it, that sheet still agrees with the physical one.
    @test amplitude(continued, 2.4 + 0.1im) ≈ amplitude(physical, 2.4 + 0.1im)
    @test amplitude(continued, 2.4) ≈ amplitude(physical, 2.4)
    lower = amplitude(continued, 2.4 - 0.1im)
    @test !(lower ≈ amplitude(physical, 2.4 - 0.1im))
    @test lower ≈ transpose(lower)
    @test isfinite(detD(continued, 2.4 - 0.1im))

    between = TMatrix(K, continue_channels(channels_cm, 1.5))
    @test detD(between, 1.5 + 0.1im) ≈ detD(physical, 1.5 + 0.1im)
    @test !(detD(between, 1.5 - 0.1im) ≈ detD(physical, 1.5 - 0.1im))
    @test isfinite(detD(between, 1.5 - 0.1im))
end

@testset "Angled Chew-Mandelstam cut" begin
    for (m1, m2) in ((1.1, 1.1), (0.14, 1.87))
        ch = TwoBodyChewMandelstamChannel(m1, m2)
        th = threshold(ch)
        first_sheet = ContinuationChannel(ch, 1)
        second_sheet = ContinuationChannel(ch, 2)
        # cut_angle is in radians. 0 is mode 1, and -π/2 is mode -90.
        down = AngledCutChannel(ch, -π / 2)
        along_axis = AngledCutChannel(ch, 0)
        flat = AngledCutChannel(ch, -π)
        shallow = AngledCutChannel(ch, -π / 6)

        @test threshold(down) == th
        @test down.cut_angle == -π / 2
        @test along_axis.cut_angle == 0
        @test_throws ArgumentError AngledCutChannel(ch, 2π)
        @test_throws ArgumentError AngledCutChannel(ch, -π - 0.1)
        @test_throws ErrorException iρ(
            AngledCutChannel(TwoBodyChewMandelstamChannel(m1, m2; L=1), -π / 2),
            th + 1,
        )

        samples = (
            th + 0.4,
            th + 0.7 + 0.2im,
            th + 0.5 * cis(-π / 12),
            th + 0.5 * cis(-π / 2),
            th + 0.5 * cis(-5π / 6),
            th - 0.4 + 0.3im,
            0.6 * th - 0.2im,
        )
        for m in samples
            @test iρ(down, m) == iρ(ContinuationChannel(ch, -90), m)
            @test iρ(along_axis, m) == iρ(ContinuationChannel(ch, 1), m)
            @test iρ(flat, m) == iρ(ContinuationChannel(ch, 12), m)
        end

        # Inside the -π/6 wedge the shallow cut is on sheet II.
        # Below that ray it is still on sheet I, unlike mode -90.
        inside = th + 0.8 * cis(-π / 12)
        below_ray = th + 0.8 * cis(-π / 2)
        @test iρ(shallow, inside) == iρ(second_sheet, inside)
        @test iρ(shallow, below_ray) == iρ(first_sheet, below_ray)
        @test iρ(shallow, th + 0.8 + 0.2im) == iρ(first_sheet, th + 0.8 + 0.2im)

        # The real-axis cut is glued; the discontinuity sits on the ray.
        m_cut = th + 1.0
        ε = 1e-6
        @test abs(iρ(shallow, m_cut + im * ε) - iρ(shallow, m_cut - im * ε)) < 1e-4
        α = -π / 6
        radius = 0.8
        toward_axis = th + radius * cis(α + 1e-3)
        past_ray = th + radius * cis(α - 1e-3)
        @test iρ(shallow, toward_axis) == iρ(second_sheet, toward_axis)
        @test iρ(shallow, past_ray) == iρ(first_sheet, past_ray)
        @test abs(iρ(shallow, toward_axis) - iρ(shallow, past_ray)) > 0.1

        # A positive angle mirrors the wedge into the upper half-plane.
        raised = AngledCutChannel(ch, π / 6)
        upper = th + 0.8 * cis(π / 12)
        lower = th + 0.8 * cis(-π / 12)
        @test iρ(raised, upper) == iρ(second_sheet, upper)
        @test iρ(raised, lower) == iρ(first_sheet, lower)
    end

    channels_cm = SVector(
        TwoBodyChewMandelstamChannel(1.5, 1.5),
        TwoBodyChewMandelstamChannel(0.5, 0.5),
        TwoBodyChewMandelstamChannel(1.0, 1.0),
    )
    # Thresholds are 3, 1, 2. The closed channel uses angle 0 so the vector
    # stays a single concrete element type.
    rotated = SVector(
        AngledCutChannel(channels_cm[1], 0),
        AngledCutChannel(channels_cm[2], -π / 6),
        AngledCutChannel(channels_cm[3], -π / 6),
    )
    @test eltype(rotated) <: AngledCutChannel
    @test getproperty.(rotated, :cut_angle) == SVector(0, -π / 6, -π / 6)
    @test iρ(rotated[1], 2.4 - 0.2im) == iρ(ContinuationChannel(channels_cm[1], 1), 2.4 - 0.2im)

    K = KMatrix([(M=3.5, gs=[1.0, 0.5, 0.2])])
    physical = TMatrix(K, channels_cm)
    continued = TMatrix(K, rotated)
    @test amplitude(continued, 2.4) ≈ amplitude(physical, 2.4)
    @test amplitude(continued, 2.4 + 0.1im) ≈ amplitude(physical, 2.4 + 0.1im)
    lower = amplitude(continued, 2.4 - 0.05im)
    @test !(lower ≈ amplitude(physical, 2.4 - 0.05im))
    @test lower ≈ transpose(lower)
    @test isfinite(detD(continued, 2.4 - 0.05im))
end
