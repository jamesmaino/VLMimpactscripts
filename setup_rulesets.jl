using FieldMetadata, GeoData, Dates

import FieldMetadata: @relimits, limits, @reflattenable, flattenable


# Override default parameter range limits
@relimits struct HumanDispersal
    human_exponent  | (0.0, 3.0)
    dist_exponent   | (0.0, 3.0)
    dispersalperpop | (0.0, 1e-8)
    max_dispersers  | (1e1, 1e4)
end
# Constant growth
@relimits struct ExactLogisticGrowth
    intrinsicrate | (0.0, 10.0)
end
# Alee
@relimits struct AlleeExtinction
    minfounders | (2e0, 1e3)
end

# Dont parametrise carrycap
@reflattenable struct ExactLogisticGrowthMap
    carrycap | false
end
# Don't parametrise λ, it will just be the top limit: we need to estimate it
@reflattenable struct ExponentialKernel
    λ | false
end

FloatType = Float32

floatconvert(a::AbstractArray) = convert(Array{FloatType,2}, a)
floatconvert(x::Number) = convert(FloatType, x)

setup_rulesets(datafile, starttime, timestep, stoptime) = begin

    # Load data
    data = h5open(datafile, "r")
    dimz = Lat((20, 60), nothing, GeoData.Ordered(GeoData.Forward(), GeoData.Reverse())), Lon((-125, -75))

    init = GeoArray(floatconvert.(read(data["x_y_initial"]["brisbane"]) .* 1e9), dimz; missingval=NaN32) # Arbitrary initial condition
    monthsofyear = Time(DateTime(2008):Month(1):DateTime(2008)+Month(11))
    popgrowth = GeoArray(replace(floatconvert.(read(data["x_y_month_intrinsicGrowthRate"])), NaN32=>0), (dimz..., monthsofyear);
                         missingval=NaN32)
    human_pop = GeoArray(replace(floatconvert.(read(data["x_y_popdens"])), NaN=>missing), dimz; missingval=missing)

    # Define simulation settings
    cellsize = floatconvert(1.0)
    framesperstep = 12
    masklayer = mask(GeoArray(popgrowth, data=(floatconvert.(read(data["x_y_month_intrinsicGrowthRate"]))))[Time(1)])


    # Rules ###########################################################

    # Human dispersal
    scale = 8
    aggregator = mean
    human_exponent = floatconvert(2.0)
    dist_exponent = floatconvert(2.0)
    dispersalperpop = floatconvert(1e-9)
    max_dispersers = floatconvert(500.0)
    shortlist_len = 100
    @time humandisp = HumanDispersal(parent(human_pop); scale=scale, shortlist_len=shortlist_len, dispersalperpop=dispersalperpop,
                                     max_dispersers=max_dispersers, cellsize=cellsize, human_exponent=human_exponent,
                                     dist_exponent=dist_exponent, timestep=timestep);

    # Climate driven growth
    carrycap = floatconvert(1e8)

    # Convert growth arrays to units
    growth = ExactLogisticGrowthMap(layer=popgrowth, carrycap=carrycap, timestep=Day(1));

    # Constant growth
    constant_growth = ExactLogisticGrowth(intrinsicrate=floatconvert(0.1), timestep=Day(1), carrycap=carrycap)

    # Local dispersal
    λ = floatconvert(0.0125)
    radius = 1
    sze = 2radius + 1

    distancemeth = AreaToArea(30)
    @time hood = DispersalKernel{radius}(;kernel=zeros(FloatType, radius, radius), cellsize=cellsize,
                                   formulation=ExponentialKernel(λ), distancemethod=distancemeth)
    localdisp = InwardsPopulationDispersal(hood)
    display(hood.kernel * carrycap)

    # Allee effects
    minfounders = floatconvert(10.0)
    allee = AlleeExtinction(minfounders=minfounders)


    # Define combinations for comparison  ##########################
    kwargs = (init=parent(init), mask=masklayer, timestep=timestep)
    full = Ruleset(humandisp, Chain(localdisp, allee, growth); kwargs...)
    nolocal = Ruleset(humandisp, allee, growth; kwargs...);
    noallee = Ruleset(humandisp, Chain(localdisp, growth); kwargs...);
    nohuman = Ruleset(Chain(localdisp, allee, growth); kwargs...);
    noclimate = Ruleset(humandisp, Chain(localdisp, allee, constant_growth); kwargs...)
    ruleset = Ruleset(Chain(localdisp, growth); kwargs...);

    ((full=full, nolocal=nolocal, noallee=noallee, nohuman=nohuman,
      noclimate=noclimate), init)
end
