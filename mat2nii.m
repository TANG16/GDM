function mat2nii(map,prefix)
%function to write matlab map to a nifti file
%Input:
% 1) map: matlab map
% 2)template: a file to borrow nifti header from, one that matches the map
% dimensions
% 3) prefix: prefix of output
%Output:
% prefix.nii
% Example usage:
%mat2nii(map{1}.stat{1},'/cbica/software/external/fsl/4.1.5/data/standard/MNI152lin_T1_2mm_brain.nii.gz','output')

template = '/sbia/sbiaprj/ADNI/ADNI_2015/Protocols/RAVENS_norm_DS94_s8/002_S_0619/2006-06-01/002_S_0619_2006-06-01_MPRAGE_scaled_LPS_brain_mars-ss_IC_QCed_MICO_C0.8-G1.0-W1.2-L10_RAVENS_150_norm_DS94_s8.nii.gz';

tmp=load_untouch_nii_gz(template);
tmp.hdr.dime.datatype=16;
tmp.img=map;

save_untouch_nii(tmp,prefix);
gzip([prefix '.nii'])
delete([prefix '.nii'])
disp(['Saved as ' prefix '.nii.gz'])
end