module Codes

using OffsetArrays 

export binarydot, poly2trellis, Trellis, convenc, vitdec 

struct Trellis
    reglen :: Int
    states :: Array{UInt8, 2}
    outputs :: Array{UInt8, 2}
    codes :: Array{UInt8, 2}
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
        z = z >> 1;
    end
    return out;
end

function foo(reglen::Integer, codes :: Vector{<:Integer})

    statelen = reglen - 1;
    nstates = 1 << statelen;

    states = OffsetArray(zeros(Integer, nstates, 2), 0:nstates-1, 0:1);
    outputs = OffsetArray(zeros(Integer, nstates, 2), 0:nstates-1, 0:1);

    mask1 = nstates;
    mask0 = 0;

    for state in 0:nstates-1
        nxt0, nxt1 = mask0 | state, mask1 | state;
        out0, out1 = 0, 0; 
        for (i, code) in enumerate(codes)
            shft = i-1;
            out0 |= (binarydot(code, nxt0) << shft);
            out1 |= (binarydot(code, nxt1) << shft);
        end

        nxt0, nxt1 = nxt0 >> 1, nxt1 >> 1; 
        outputs[state, :] = [out0, out1];
        states[state, :] = [nxt0, nxt1];
    end

    return parent(outputs), parent(states);
end

# function convenc(trellis::Trellis, bits::Array{UInt8, 2}) :: Array{UInt8,2}
#     state = 0x00;
#     output = zeros(UInt8, 1, length(bits));
#     for (i, bit) in enumerate(bits)
#         row, col = state + 1, bit+1;
#         output[i] = trellis.outputs[row, col];
#         state = trellis.states[row, col];
#     end
#     return output;
# end

# function minsum(errors, observed, predictions) 
#     n = length(errors);
#     half = n >> 1 

#     e = zeros(1, n);
#     p = zeros(1, n);

#     for i in 0:n-1
#         if i <= half
#             i0 = ((i-1) << 1) + 1
#             i1 = i0 + 1;
#             e0 = errors[i0] + biterror(observed, predictions[i0, 1]);
#             e1 = errors[i1] + biterror(observed, predictions[i1, 1]);
#         else
#             i0 = ((i - 1 - half) << 1) + 1;
#             i1 = i0 + 1;
#             e0 = errors[i0] + biterror(observed, predictions[i0, 2]);
#             e1 = errors[i1] + biterror(observed, predictions[i1, 2]);            e0 = errors[i0+1] + 
#         end
        
#         if e0 < e1 
#             p[i], e[i] = i0, e0;
#         else 
#             p[i], e[i] = i1, e1;
#         end
#     end

#     return e, p
# end

end # module

