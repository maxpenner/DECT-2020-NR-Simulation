function [] = check_equality_of_two_numbers(a,b,error_margin)

    if error_margin == 0
        if a ~= b
            error('Numbers are not exactly equal: a=%f b=%f', a, b);
        end
    else
        if a>b
            c = b/a;
        else
            c = a/b;
        end
        if abs(c-1) > error_margin
            error('Difference not within the error margin: a=%f b=%f error_margin=%f', a, b, error_margin);
        end
    end
end

