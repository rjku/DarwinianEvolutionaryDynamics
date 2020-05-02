
module mPlot

using PyPlot

function draw_boxplot(data::Vector{<:AbstractVector},bins::AbstractVector,width::Real,offset::Real,fill_color="tab:blue";showOutliers::Bool=false)
    pst = bins[2:end] .+ offset .- (bins[2] - bins[1]) / 2.
    bp = boxplot(data, positions=pst, widths=width, patch_artist=true, showfliers=showOutliers) # manage_ticks=false, 
    for patch in bp["boxes"]
        patch.set(facecolor=fill_color)
    end
    for patch in bp["medians"]
        patch.set(color="black")
    end
end

export draw_boxplot

end  # module mPlot