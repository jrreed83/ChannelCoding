module Codes

using OffsetArrays 

export binarydot, poly2trellis, Trellis, convenc, vitdec, minsum

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

    state2idx = Dict([(i-1:i) for i in 1:nstates]);

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

        idx = state2idx[state];

        outputs[idx, :] = [out0, out1];
        states[idx, :] = [nxt0, nxt1];
    end

    return Trellis(nstates, states, outputs);
end

#
function convenc(t::Trellis, x::Vector{<:Integer}) :: Vector{Integer}
    y = zeros(Integer, length(x));
    n = t.nstates;

    state2idx = Dict([(i-1:i) for i in 1:n]);    
    bit2idx = Dict(0=>1, 1=>2);
    h = 0;
    for (i, xi) in enumerate(x)
        row = state2idx[h];
        col = bit2idx[xi]; 
        y[i] = t.outputs[row, col];
        h = t.states[row, col];
    end
    return y;
end

function minsum(t :: Trellis, errors :: Vector{Float64}, observed :: Integer) 
    n = length(errors);
    half = n >> 1 

    e = zeros(Float64, n);
    p = zeros(Integer, n);

    preds = t.outputs;

    state2idx = Dict([(i-1:i) for i in 1:n]);    
    bit2idx = Dict(0=>1, 1=>2);

    for i in 0:n-1
        if i < half
            i0 = i << 1;
            i1 = i0 + 1;
            e0 = errors[state2idx[i0]] + biterror(observed, preds[state2idx[i0], bit2idx[0]]);
            e1 = errors[state2idx[i1]] + biterror(observed, preds[state2idx[i1], bit2idx[0]]);
        else
            i0 = (i-half) << 1;
            i1 = i0 + 1;
            e0 = errors[state2idx[i0]] + biterror(observed, preds[state2idx[i0], bit2idx[1]]);
            e1 = errors[state2idx[i1]] + biterror(observed, preds[state2idx[i1], bit2idx[1]]);            
        end
        
        p[state2idx[i]], e[state2idx[i]] = e0 < e1 ? (i0, e0) : (i1, e1);
    end

    return e, p;
end

function vitdec(t::Trellis, inputs::Vector{<:Integer})

    n = length(inputs);
    paths = zeros(Integer, t.nstates, n);

    state2idx = Dict([(i-1,i) for i in 1:n]);  
    idx2state = Dict([(i,i-1) for i in 1:n]);    

    errs = fill(Inf, t.nstates);
    errs[state2idx[0]] = 0;
    
    # Forward pass
    for (i, obs) in enumerate(inputs)
        errs, paths[:,i] = minsum(t, errs, obs);
    end

    # Backward pass 
    bits = zeros(Integer, n);    
    k1 = idx2state[argmin(errs)];
    for i in n:-1:1
        k0 = paths[state2idx[k1], i];
        bits[i] = (k1 == (k0 >> 1)) ? 0 : 1;
        k1 = k0;
    end

    return bits;
end

end # module

