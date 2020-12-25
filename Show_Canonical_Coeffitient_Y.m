function [] = Show_Canonical_Coeffitient_Y(CC , name , Data_DIR)

MNI_1 = 60;
MNI_2 = 72;
MNI_3 = 60;

load(Data_DIR);

mask = zeros(MNI_1 , MNI_2 , MNI_3);
for k=1:MNI_3
    for i=1:MNI_1
        for j=1:MNI_2
            if (Tagged_Obj_Y(i,j,k) ~= 0)
                mask(i,j,k) = CC(Tagged_Obj_Y(i,j,k),1);
            end
        end
    end
end

% mask = zeros(MNI_1 , MNI_2 , MNI_3);
% 
% out_data = reshape(CC , D1 , D2 , D3);
% out_data = permute(out_data, [2 , 1 , 3]);
% 
% %reshape to MNI 2mm space 91*109*91
% mask(Min_X : Min_X + D2 - 1 , Min_Y : Min_Y + D1 - 1 , Min_Z : Min_Z + D3 - 1) = out_data(: , : , :);  

out_data_nii = make_nii(mask);

out_name = strcat(pwd,name);
save_nii(out_data_nii , out_name);

new_data = load_nii(out_name);
view_nii(new_data);

