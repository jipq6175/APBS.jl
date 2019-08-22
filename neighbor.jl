function neighbor(photonstream::Array{Int64}, lb::Int64, ub::Int64)
    g::Int64 = length(filter(x-> x>lb && x<ub, photonstream));
    return g
end
