using CairoMakie
using ScatteringKMatrix
using StaticArrays

const mπ = 0.14
const resonance_mass = 0.77
const target_width = 0.15

# For T = g²/(M²-m²-g² iρ), the narrow-width estimate is
# Γ ≈ g²ρ(M)/M. This gives a width close to 0.15 for the square-root
# model and provides one common coupling for comparing the two phase spaces.
const ρ_on_shell = imag(iρ(TwoBodyChannel(mπ, mπ), resonance_mass))
const gsq = target_width * resonance_mass / ρ_on_shell
const K = KMatrix([(M=resonance_mass, gs=[sqrt(gsq)])])

const channels = (
    ("Chew–Mandelstam", TwoBodyChewMandelstamChannel(mπ, mπ)),
    ("Square root", TwoBodyChannel(mπ, mπ)),
)
const modes = (("sheet I", 1), ("sheet II", 2), ("downward cut", -90))

const re_m = range(0.12, 1.12; length=401)
const im_m = range(-0.50, 0.50; length=401)

function imaginary_amplitude(channel, mode)
    continued = SVector(ContinuationChannel(channel, mode))
    amplitude_model = TMatrix(K, continued)
    return [
        imag(amplitude(amplitude_model, complex(x, y))[1, 1])
        for x in re_m, y in im_m
    ]
end

values = [imaginary_amplitude(channel, mode) for (_, channel) in channels, (_, mode) in modes]
finite_values = filter(isfinite, reduce(vcat, vec.(values)))
sorted_magnitudes = sort!(abs.(finite_values))
color_limit = sorted_magnitudes[round(Int, 0.985 * length(sorted_magnitudes))]

set_theme!(Theme(
    fontsize=19,
    Axis=(
        backgroundcolor=:transparent,
        xgridvisible=false,
        ygridvisible=false,
    ),
))

figure = Figure(size=(1500, 940), backgroundcolor=:white)
plots = Matrix{Any}(undef, length(channels), length(modes))

for (row, (channel_name, _)) in enumerate(channels)
    for (column, (mode_name, mode)) in enumerate(modes)
        axis = Axis(
            figure[row, column];
            title=row == 1 ? "$(mode_name)  (mode $(mode))" : "",
            xlabel=row == length(channels) ? "Re m" : "",
            ylabel=column == 1 ? "$(channel_name)\nIm m" : "",
            aspect=DataAspect(),
        )
        plots[row, column] = heatmap!(
            axis,
            re_m,
            im_m,
            values[row, column];
            colormap=:vik,
            colorrange=(-color_limit, color_limit),
            rasterize=true,
        )

        threshold_mass = 2mπ
        if mode in (1, 2)
            lines!(axis, [threshold_mass, last(re_m)], [0, 0]; color=(:black, 0.65), linestyle=:dash, linewidth=2)
        else
            lines!(axis, [threshold_mass, threshold_mass], [0, first(im_m)]; color=(:black, 0.65), linestyle=:dash, linewidth=2)
        end
        scatter!(axis, [threshold_mass], [0]; color=:black, markersize=8)
        scatter!(axis, [resonance_mass], [0]; color=:white, strokecolor=:black, strokewidth=1.5, markersize=10)
        xlims!(axis, extrema(re_m))
        ylims!(axis, extrema(im_m))
    end
end

Colorbar(
    figure[:, end + 1],
    plots[1, 1];
    label="Im T(m)",
    width=24,
)
Label(
    figure[0, :],
    "Single-channel continuation: m₁=m₂=$(mπ), M=$(resonance_mass), " *
    "g²=$(round(gsq; digits=3)) (Γ≈$(target_width))";
    fontsize=24,
    font=:bold,
)

output = joinpath(@__DIR__, "continuation_heatmaps.png")
save(output, figure; px_per_unit=1.5)
println("saved $(output); g²=$(gsq), color limit=$(color_limit)")
