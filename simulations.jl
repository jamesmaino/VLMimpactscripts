# Ativate the Julia project

# Load dependencies
# using Pkg
# Pkg.activate(".")
using Revise, GeoData, ArchGDAL; const AG = ArchGDAL
using GeoData: Lat, Lon, Time
using HDF5, JLD2, Dates, DynamicGrids, Dispersal, CSV, Statistics
using Plots, Distributions


# Load simualtions for USA
datafile = "spread_inputs_Aus_SWD.h5"
starttime = DateTime(2008, 6)
timestep = Month(1)
stoptime = DateTime(2018, 12)
tspan = starttime, stoptime
include("setup_rulesets.jl")
sim_rulesets, init = setup_rulesets(datafile, starttime, timestep, stoptime);

ruleset = sim_rulesets[:full];

output = ArrayOutput(init, length(starttime:timestep:stoptime))
sim!(output, ruleset; init=init, tspan=tspan)

cropval = AG.registerdrivers() do
    AG.read("data/clum_crop_value_cea.tif") do dataset
    band1 = AG.getband(dataset,1)
    new1 = AG.read(band1)
    end
end;
cropval = transpose(cropval)
cropval = map(x -> x < 0 ? 0.0 : x, cropval)
sum(cropval)
hort =  map(x -> x > 0 ? 1 : 0, cropval)

heatmap(reverse(log.(cropval .+ 1), dims = 1))
heatmap(reverse(convert(Array, output[100]), dims = 1))

# build trap network
trappinggrid = similar(cropval, Int8);
fill!(trappinggrid, 0);
traps_per_cell = 0.1 # how many grid cells between traps
for i = 1:length(trappinggrid)
    if cropval[i] > 0
        trappinggrid[i] =
        rand(Poisson(traps_per_cell))
    end
end
sum(trappinggrid)
heatmap(reverse(trappinggrid, dims=1))

p_detection(detection_rate_m, popdens, ntraps; trap_coverage = 2.7) =
    1 - ((1 - detection_rate_m)^ntraps)^(popdens * trap_coverage)

# findall(x -> x >0, output[68])

# get suitability
data = h5open(datafile, "r")
popgrowth = replace(floatconvert.(read(data["x_y_month_intrinsicGrowthRate"])), NaN32=>0);
popgrowth = [mean(popgrowth[i,j,:]) for i = 1:size(popgrowth)[1], j = 1:size(popgrowth)[2]];
suit = map(x -> x > 0 ? 1 : 0, popgrowth)
heatmap(reverse(suit, dims = 1))
suithort = suit .* hort

function area_at_detection(output, trappinggrid)
    gridsize = 9 * 9 * 100 # ha per 9 km gridcell
    detection_rate = 0.0052 # mean rate of detection for 93 m trapping radius from Kirkpatrick 2018 across mean([27, 48]) = 37.5 d
    detection_rate_m = detection_rate/37.5*30
    for timestep = 1:length(output)
        for cell = 1:length(output[timestep])
            popdens = output[timestep][cell]/gridsize # flies per ha
            trapnumber = trappinggrid[cell]
            if trapnumber == 0 || popdens == 0
                continue
            end
            p = p_detection(detection_rate_m, popdens, trapnumber, trap_coverage = 2.7)
            if rand(Binomial(1, p)) == 1
                return (timestep, sum(output[timestep] .> 0) * gridsize)
            end
        end
    end
    return missing # not detected
end

function build_trappinggrid(traps_per_cell)
    trappinggrid = similar(cropval, Int8)
    fill!(trappinggrid, 0)
    for i = 1:length(trappinggrid)
        if cropval[i] > 0
            trappinggrid[i] =
            rand(Poisson(traps_per_cell))
        end
    end
    return trappinggrid
end

spread_reps = 10
traps_per_cells = 1 ./ 10 .^ (-1:0.1:1)
surveillance_reps = 10

time_area = [(-1, -1) for i = 1:spread_reps, j = 1:length(traps_per_cells), k = 1:surveillance_reps]
for i = 1:spread_reps
    # initialise in suitable horticulture cell
    init = reshape(rand(Multinomial(1, vec(suithort ./ sum(suithort)))), size(hort)) * 1e9
    output = ArrayOutput(init, length(starttime:timestep:stoptime))
    sim!(output, ruleset; init=init, tspan=tspan)
    for j = 1:length(traps_per_cells)
        for k = 1:surveillance_reps
            println([i, j, k])
            trappinggrid = build_trappinggrid(traps_per_cells[j])
            time_area[i, j, k] = area_at_detection(output, trappinggrid)
            println(time_area[i, j, k] )
        end
    end
end

using DataFrames
d = DataFrame([Float32,Float32,Float32,Float32],[:time, :timesd, :area, :areasd], length(traps_per_cells))
d.traps_per_cell = traps_per_cells

for j = 1:length(traps_per_cells)
    println()
    area_invaded = map(x->x[2], time_area[:,j,:])
    detection_time = map(x->x[1], time_area[:,j,:])
    d[j, :area] = mean(skipmissing(area_invaded))
    d[j, :areasd] = std(skipmissing(area_invaded))
    d[j, :time] = mean(skipmissing(detection_time))
    d[j, :timesd] = std(skipmissing(detection_time))
end

plot(d.traps_per_cell, (d.area), yerror =d.area, xlabel = "trap per cell", label = "area invaded at detection")

plot(d.traps_per_cell, log10.(d.area), xlabel = "trap per cell", label = "area invaded at detection")
