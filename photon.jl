# Read photonstream and convertphotonstream in julia
# bitbyte: 8bit or 16bit
# systemclock: 40 or 60 MHz
# arrivaltime = false: only one output of clockticks
function photon(file::String, arrivaltime::Bool = false)

    # Read in the headers
    f = open(file, "r");
    bitbyte = read(f, UInt8);
    systemclock = read(f, UInt8);

    # Identify the data format
    if bitbyte == 8
        T = UInt8;
        warn("Uisng 8-bit data format.");
    elseif bitbyte == 16
        T = UInt16;
        #info("Using 16-bit data format.");
    else
        error("Error: The data format is not specified.");
    end

    # Read in the binary file
    rawdata = Array{T}(0);
    while true
        try
            n = read(f, T);
            push!(rawdata, n);
        catch errmsg
            if isa(errmsg, EOFError);
                break;
            else
                throw(errmsg);
            end
        end
    end
    flush(f);
    close(f);

    # Data processing based on the format
    # This part is tricky: case 1: [65535 65535 5598]
    #                      case 2: [65535 5598 6633]
    if bitbyte == 8
        warn("Still working on this....");
    elseif bitbyte == 16
        data = rawdata + 0;
        index = find(indexin(rawdata, [65535]));
        consecutive = find(indexin(diff(index), [1]));
        cindex = Array{Int64}(0);

        #= Deal with consecutive cases (Case 1)
        for i in consecutive
            data[index[i]:(index[i]+2)] = 0;
            c = bin(rawdata[index[i]+2], 16) * bin(rawdata[index[i]+1], 16);
            data[index[i]+1] = parse(BigInt, c, 2);
            push!(cindex, index[i], index[i]+1);
        end
        =#

        # Still deal with consecutive case but this is easier
        for i in consecutive
            push!(cindex, index[i]+1);
        end

        # Deal with the full words (Case 2)
        filter!(x->!(x in cindex), index);
        for i in index
            c = bin(rawdata[i+2], 16) * bin(rawdata[i+1], 16);
            data[i] = parse(BigInt, c, 2);
            data[i+1:i+2] = 0;
        end
        filter!(x->x≠0, data);
    end

    # Convert photon data to arrival time
    timebin = 1.0/(systemclock * 1e6);
    maxphoton = 2^bitbyte - 1;
    id = indexin(data, [maxphoton]);
    i = ~convert(BitArray, id);
    cummtime = i .* cumsum(data);
    filter!(x -> x!=0, cummtime);
    timedata = cummtime * timebin;

    # Determine if arrivaltime is in output
    if arrivaltime
        return data, timedata;
    else
        return data;
    end
end
