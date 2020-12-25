function Write_CV_data_2_Excel(U , V , training_sample , N_G1 , N_G2 , N_train_G1 , i_CV , fname , sname_1 , sname_2 , ...
            range_1 , range_2 , range_3 , range_4)

%     data_1 = G(training_sample(U(1 : N_train_G1 , 1) <  0));       %The HC subjects that their CV is negative
%     data_2 = setdiff(G(training_sample(U( : , 1) >  0)) , G_1);      %The AD subjects that their CV is positive
data_1 = training_sample(U(1 : N_train_G1 , 1) <  0);       %The HC subjects that their CV is negative
data_2 = setdiff(training_sample(U( : , 1) >  0) , training_sample(1 : N_train_G1)) - N_G1;
    
%     data_3 = G(training_sample(V(1 : N_train_G1 , 1) <  0));       %The HC subjects that their CV is negative
%     data_4 = setdiff(G(training_sample(V( : , 1) >  0)) , G_1);      %The AD subjects that their CV is positive
data_3 = training_sample(V(1 : N_train_G1 , 1) <  0);       %The HC subjects that their CV is negative
data_4 = setdiff(training_sample(V( : , 1) >  0) , training_sample(1 : N_train_G1)) - N_G1;

temp_1 = zeros(1 , N_G1);
temp_2 = zeros(1 , N_G2);
temp_3 = zeros(1 , N_G1);
temp_4 = zeros(1 , N_G2);

temp_1(data_1) = 1;
temp_2(data_2) = 1;
temp_3(data_3) = 1;
temp_4(data_4) = 1;

if i_CV == 1
    xlswrite(fname , temp_1 , sname_1 , range_1);  % write the new data after the old data
    xlswrite(fname , temp_2 , sname_2 , range_2);  % write the new data after the old data
    xlswrite(fname , temp_3 , sname_1 , range_3);  % write the new data after the old data
    xlswrite(fname , temp_4 , sname_2 , range_4);  % write the new data after the old data
else
    num = xlsread(fname , sname_1 , range_1);
    xlswrite(fname , temp_1 + num , sname_1 , range_1);  % write the new data after the old data

    num = xlsread(fname , sname_2 , range_2);
    xlswrite(fname , temp_2 + num , sname_2 , range_2);  % write the new data after the old data

    num = xlsread(fname , sname_1 , range_3);
    xlswrite(fname , temp_3 + num , sname_1 , range_3);  % write the new data after the old data

    num = xlsread(fname , sname_2 , range_4);
    xlswrite(fname , temp_4 + num , sname_2 , range_4);  % write the new data after the old data
end
        