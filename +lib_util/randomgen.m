function out_d = randomgen(reset_internal_state)

    persistent x;
    persistent y;
    persistent z;
    persistent w;
    
    if nargin < 1
        reset_internal_state = false;
    end

    if isempty(x) || reset_internal_state==true
        x=uint32(0);%123456789);
        y=uint32(362436069);
        z=uint32(521288629);
        w=uint32(88675123);
    end
    
    if reset_internal_state==true
        x=uint32(0);%123456789);
        y=uint32(362436069);
        z=uint32(521288629);
        w=uint32(88675123);
        
        out_d = 0;
        return;
    end

    t = uint32(bitxor(x, bitshift(x,11)));
    x = y;
    y = z;
    z = w;
    
    a = bitsra(w,19);
    b = bitsra(t,8);
    c = bitxor(t,b);
    
    w = bitxor(bitxor(w,a),c);
    
    int64_of_UINT32_MAX_half = bitsra(int64(4294967295), 1);
    
    % expand to 64 bit
    out = int64(w);
    
    % center
    out = out - int64_of_UINT32_MAX_half;
    
    % normalize to +/- 1
    out_d = double(out)/double(int64_of_UINT32_MAX_half);
end