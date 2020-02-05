
module mPlot

using PyPlot

function draw_boxplot(data::Vector{Vector{Float64}},bins::AbstractVector,width::Real,offset::Real,fill_color="tab:blue";showOutliers::Bool=false)
    pst = bins[2:end] .+ offset
    bp = boxplot(data, positions=pst, widths=width, patch_artist=true, manage_ticks=false, showfliers=showOutliers)
    for patch in bp["boxes"]
        patch.set(facecolor=fill_color)
    end
    for patch in bp["medians"]
        patch.set(color="black")
    end
end

export draw_boxplot

end  # module mPlot
