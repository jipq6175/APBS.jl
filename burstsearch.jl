# FRET Burst search in Julia
# [M, L, W] = [25, 30, 0.001];
# options = "allphoton", "donor" or "acceptor"
# The performance is optimized by the integer computation
# The performance is improved by 20 times by avoiding global Variables
# The performance is 6x faster than MATLAB

@everywhere include("photon.jl");
@everywhere include("neighbor.jl");

function burstsearch(donorfile::String, acceptorfile::String, M::Int64=25, L::Int64=30, W::Float64=0.001, options::String="allphoton")

    # Fetch data from "photon.jl"
    d::Array{Int64} = photon(donorfile);
    a::Array{Int64} = photon(acceptorfile);
    allphoton = Dict("donor" => cumsum(d), "acceptor" => cumsum(a));
    allphoton["all"] = sort(vcat(cumsum(d), cumsum(a)));


    # Determine which channel to search
    if options == "allphoton"
        photonstream = allphoton["all"];
    elseif options == "donor"
        photonstream = allphoton["donor"];
    elseif options == "acceptor"
        photonstream = allphoton["acceptor"];
    else
        photonstream = allphoton["all"];
        warn("Specified option is not found. Using allphoton channel by default.");
    end

    # Burst searching algorithm
    burst = Array{Int64}(0);
    burstevents = Array{Int64}(0);
    fretphotons = Array{Int64}(0);
    et = Array{Float64}(0);

    systemclock::Int64 = 40;
    window::Int64 = W * systemclock * 1000000;
    l::Int64 = 0;
    count::Int16 = 0;
    for i in collect(1:length(photonstream))
        if photonstream[i] < window/2
            continue;
        else
            center = photonstream[i];
            lb::Int64 = center - window/2;
            ub::Int64 = center + window/2;
            m = neighbor(photonstream, lb, ub);
            if m > M
                l += 1;
            else
                if l < L
                    l = 0;
                else
                    burst = vcat(burst, photonstream[i-1-l:i]);
                    count += 1;
                    println("The $count burst has $l photons.");
                    burstevents = vcat(burstevents, [count l]);
                    # Calculate FRET
                    if options == "allphoton"
                        na = neighbor(allphoton["acceptor"], photonstream[i-1-l], photonstream[i]);
                        fretphotons = vcat(fretphotons, [l - na na]);
                        et = vcat(et, [photonstream[i-1-l]/(systemclock * 1e6) na/l]);
                        l = 0;
                    end
                end
            end
        end
    end

    # Outputs
    println("The burst searching is complete.");
    output = Dict("burstevents" => [burstevents fretphotons], "et" => et);
    return output;
end
