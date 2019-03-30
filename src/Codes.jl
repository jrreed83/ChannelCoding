module Codes

using OffsetArrays 

export binarydot, poly2trellis, Trellis, convenc, vitdec 

struct Trellis
    nstates :: Integer 
    states  :: Array{Integer, 2}
    outputs :: Array{Integer, 2}
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

    s = zeros(Integer, nstates, 2);
    o = zeros(Integer, nstates, 2);
    
    states = OffsetArray(s, 0:nstates-1, 0:1);
    outputs = OffsetArray(o, 0:nstates-1, 0:1);

    mask1 = nstates;
    mask0 = 0;

    ncodes = length(codes); 

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

        outputs[state, :] = [out0, out1];
        states[state, :] = [nxt0, nxt1];
    end

    return Trellis(nstates, s, o);
end

#
function convenc(trellis::Trellis, x::Vector{<:Integer}) :: Vector{Integer}
    y = zeros(Integer, length(x));
    
    n, _ = size(trellis.states);

    states = OffsetArray(trellis.states, 0:n-1, 0:1);
    outputs = OffsetArray(trellis.outputs, 0:n-1, 0:1);

    h = 0;
    for (i, xi) in enumerate(x)
         y[i] = outputs[h, xi];
         h = states[h, xi];
     end
     return y;
end

function minsum(t :: Trellis, errors :: Vector{<:Integer}, observed :: Integer) 
    n = length(errors);
    half = n >> 1 

    e = zeros(Integer, n);
    p = zeros(Integer, n);

    predictions = OffsetArray(trellis.outputs, 0:n-1, 0:1);
    errors = OffsetArray(errors, 0:n-1);
    
    for i in 0:n-1
        if i < half
            i0 = i << 1;
            i1 = i + 1;
            e0 = errors[i0] + biterror(observed, predictions[i0, 0]);
            e1 = errors[i1] + biterror(observed, predictions[i1, 0]);
         else
            i0 = (i-half) << 1
            i1 = i0 + 1;
            e0 = errors[i0] + biterror(observed, predictions[i0, 1]);
            e1 = errors[i1] + biterror(observed, predictions[i1, 1]);            
         end
        
#         if e0 < e1 
#             p[i], e[i] = i0, e0;
#         else 
#             p[i], e[i] = i1, e1;
#         end
#     end

#     return e, p
# end

end # module

