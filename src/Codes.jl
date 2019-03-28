module Codes

export binarydot, poly2trellis 

function binarydot(x::Unsigned, y::Unsigned) :: UInt8
    z = Unsigned(x) & Unsigned(y);
    out = 0x00;  
    while z != 0
        out = xor(out, z & 1);
        z >>= 1
    end
    return UInt8(out);
end

function poly2trellis(reglen::Integer, codes :: Array{UInt8,2})

    statelen = reglen - 1;
    nstates = 1 << statelen;

    nxtstates = zeros(UInt8, nstates, 2);
    outputs = zeros(UInt8, nstates, 2);

    mask1 = UInt8(nstates);
    mask0 = 0x00;


    for idx in range(1, nstates, step=1)
        state = UInt8(idx-1);
        nxt0, nxt1 = mask0 | state, mask1 | state
        out0, out1 = 0x00, 0x00; 
        for (i, code) in enumerate(codes)
            shft = i-1;
            out0 |= (binarydot(code, nxt0) << shft);
            out1 |= (binarydot(code, nxt1) << shft);
        end
        nxt0, nxt1 = nxt0 >> 1, nxt1 >> 1; 
        outputs[idx, :] = [out0, out1];
        nxtstates[idx, :] = [nxt0, nxt1];
    end

    return nxtstates, outputs
end

end # module
