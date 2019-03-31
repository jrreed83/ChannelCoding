module Codes


export 
    binarydot, 
    poly2trellis, 
    Trellis, 
    convenc, 
    vitdec, 
    minsum, 
    biterror 


state2idx(n) = Dict(i-1 => i for i in 1:n);  
idx2state(n) = Dict(i => i-1 for i in 1:n);
bit2idx = Dict(0=>1, 1=>2);

struct Trellis
    nstates   :: Integer 
    states    :: Array{Integer, 2}
    outputs   :: Array{Integer, 2}
end

function binarydot(x::Integer, y::Integer) :: Integer
    z = x & y;
    out = 0;  
    while z != 0
        out = xor(out, z & 1);
        z >>= 1
    end
    return out;
end

function biterror(x::Integer, y::Integer) :: Integer 
    z = xor(x, y);
    out = 0;
    while z != 0
        out += (z & 1);
        z >>= 1;
    end
    return out;
end

function poly2trellis(reglen::Integer, codes :: Vector{<:Integer})

    statelen = reglen - 1;
    nstates = 1 << statelen;

    states = zeros(Integer, nstates, 2);
    outputs = zeros(Integer, nstates, 2);

    mask1 = nstates;
    mask0 = 0;

    ncodes = length(codes); 

    s2i = state2idx(nstates);

    for state in 0:nstates-1
        nxt0, nxt1 = mask0 | state, mask1 | state;
        out0, out1 = 0, 0; 
        for (i, code) in enumerate(codes)
            shft = i - 1;
            out0 |= (binarydot(code, nxt0) << shft);
            out1 |= (binarydot(code, nxt1) << shft);
        end

        # Shift
        nxt0 >>= 1;
        nxt1 >>= 1;

        outputs[s2i[state], :] = [out0, out1];
        states[s2i[state], :] = [nxt0, nxt1];
    end

    return Trellis(nstates, states, outputs);
end

#
function convenc(t::Trellis, x::Vector{<:Integer}) :: Vector{Integer}
    y = zeros(Integer, length(x));
    n = t.nstates;

    s2i = state2idx(n);    
    b2i = bit2idx;
    h = 0;
    for (i, bit) in enumerate(x) 
        y[i] = t.outputs[s2i[h], b2i[bit]];
        h = t.states[s2i[h], b2i[bit]];
    end
    return y;
end

function minsum(t :: Trellis, errors :: Vector{Float64}, observed :: Integer) 
    n = length(errors);
    half = n >> 1 

    e = zeros(Float64, n);
    p = zeros(Integer, n);

    preds = t.outputs;

    s2i = state2idx(n);    
    b2i = bit2idx;

    for i in 0:n-1
        if i < half
            i0 = i << 1;
            i1 = i0 + 1;
            e0 = errors[s2i[i0]] + biterror(observed, preds[s2i[i0], b2i[0]]);
            e1 = errors[s2i[i1]] + biterror(observed, preds[s2i[i1], b2i[0]]);
        else
            i0 = (i-half) << 1;
            i1 = i0 + 1;
            e0 = errors[s2i[i0]] + biterror(observed, preds[s2i[i0], b2i[1]]);
            e1 = errors[s2i[i1]] + biterror(observed, preds[s2i[i1], b2i[1]]);            
        end
        
        p[s2i[i]], e[s2i[i]] = e0 < e1 ? (i0, e0) : (i1, e1);
    end

    return e, p;
end

function vitdec(t::Trellis, inputs::Vector{<:Integer})

    n = length(inputs);
    paths = zeros(Integer, t.nstates, n);

    s2i = state2idx(t.nstates);
    i2s = idx2state(t.nstates);
    errs = fill(Inf, t.nstates);
    errs[s2i[0]] = 0;
    
    # Forward pass
    for (i, obs) in enumerate(inputs)
        errs, paths[:,i] = minsum(t, errs, obs);
    end

    # Backward pass 
    bits = zeros(Integer, n);    
    k1 = i2s[argmin(errs)];
    for i in n:-1:1
        k0 = paths[s2i[k1], i];
        bits[i] = (k1 == (k0 >> 1)) ? 0 : 1;
        k1 = k0;
    end
    return bits;
end

end # module

