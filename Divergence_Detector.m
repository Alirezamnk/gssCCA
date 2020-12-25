function [min_input , count_input , exit] = Divergence_Detector(Input , min_input , count_input , count_max)

exit = 0;
if Input < min_input
    min_input = Input;
else
    count_input = count_input + 1;
    if count_input >= count_max
        warning('It should be stoped');
        exit = 1;
        return;
    end
end
