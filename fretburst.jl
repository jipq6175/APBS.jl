# Fretburst parallel compute all the ".dat" file in a directory

function fretburst(dir::String, gamma::Float64=1.0, beta::Float64=0.0, M::Int64=25, L::Int64=30, W::Float64=0.001, options::String="allphoton")
    #curr = pwd();

    # Introduce parallel workers
    # 4 workers seem to work the best
    nw::Int64 = 4;
    num_workers = length(workers());
    if num_workers < nw
        addprocs(nw - num_workers);
        info("# of workers: from $num_workers increased to $nw.");
    elseif num_workers > nw
        p = workers();
        rmprocs(p[nw:end]);
        warn("# of workers: from $num_workers reduced to $nw.");
    else
        info("# of workers: Already exists $nw workers.");
    end

    # include burstsearch at every worker
    @everywhere include("C:\\Users\\Yen-Lin\\Google Drive\\13. Julia\\03. smFRET\\burstsearch.jl");

    cd(dir);
    dirlist::Array{String} = readdir();
    dfiles::Array{String} = sort(filter(x -> x[end-4:end] == "A.dat", dirlist));
    afiles::Array{String} = sort(filter(x -> x[end-4:end] == "B.dat", dirlist));

    # Check the agreement on the number of files
    if length(dfiles) == length(afiles)
        l::Int64 = length(dfiles);
    else
        error("The numbers of donor/acceptor files don't match.");
    end

    # Run burstsearch function
    stage = Array{Float64}(0);
    stage = @parallel vcat for i = 1:l
        cd(dir);
        searchresult = burstsearch(dfiles[i], afiles[i], M, L, W, options);
        println("The $i th data:  $(length(searchresult["et"][:,2])) bursts were found.");
        searchresult;
    end

    # Remove additional parallel workers
    rmprocs(workers());
    info("Remove all parallel workers.");

    # Rearange all the shared arrays
    bursts = Array{Int64}(0);
    et = Array{Float64}(0);

    for i = 1:l
        bursts = vcat(bursts, stage[i]["burstevents"]);
        e = stage[i]["et"];
        # FRET time trace
        e[:,1] = e[:,1] + 30*(i-1);
        et = vcat(et, e);
    end

    l = length(et[:,2]);
    println("FRET burst has found $l events.");
    cd(curr);
    output = Dict("bursts" => bursts, "et" => et);
    return output;
end
