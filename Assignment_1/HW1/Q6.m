%% CS754-2022, Assignment 1
% Arpon Basu and Shashwat Garg
% 200050013 & 200050130

% We start with adding the code for the reading the cars
%% a) We extract first T frames, and convert to grayscale
% We obtain the T frames in the first step directly, rather
% than going for T=3, 5 and 7 separately.
file_name='cars.avi';
addpath('./MMread');
T = 20;

video_data=mmread(file_name,1:T,[],false,true);

H = 120;
W = 240;

Initial_H=video_data.height;
Initial_W=video_data.width;

frame_list=zeros(H,W,T);

for index=1:T   
    frame=rgb2gray(video_data.frames(index).cdata);    
    frame_list(:,:,index)=frame(Initial_H-H:Initial_H-1,Initial_W-W:Initial_W-1);
end
% By now, we have converted the first 7 frames to grayscale

frame_list_3 = frame_list(:,:,1:3);
frame_list_5 = frame_list(:,:,1:5);
frame_list_7 = frame_list(:,:,1:7);


%% b) Generation of Random Code matrix.

Random_C=randi([0 1],H,W,T);

Random_C_3 = Random_C(:,:,1:3);
Random_C_5 = Random_C(:,:,1:5);
Random_C_7 = Random_C(:,:,1:7);

gaussian_noise_stddev = 2;
gaussian_noise_mean = 0;

E_3=zeros(H,W);
E_5=zeros(H,W);
E_7=zeros(H,W);

gaussian_noise=gaussian_noise_mean + rand([H,W])*gaussian_noise_stddev;

for i=1:3
    E_3=E_3+frame_list_3(:,:,i).*Random_C_3(:,:,i);
end
E_3 = E_3 + gaussian_noise; 
figure('name','Coded snapshot with noise for T=3');
imshow(E_3/(3*255));
title('T=3');


for i=1:5
    E_5=E_5+frame_list_5(:,:,i).*Random_C_5(:,:,i);
end
E_5 = E_5 + gaussian_noise; 
figure('name','Coded snapshot with noise for T=5');
imshow(E_5/(5*255));
title('T=5');

for i=1:7
    E_7=E_7+frame_list_7(:,:,i).*Random_C_7(:,:,i);
end
E_7 = E_7 + gaussian_noise; 
figure('name','Coded snapshot with noise for T=7');
imshow(E_7/(7*255));
title('T=7');


close all

















%% Helper functions.

function reconstructed_frames = reconstruction_omp(encoding_code, encoded_frames, size_patch, error_limit)
    dct_2d_psi_matrix = kron(dcrmtx(size_patch)',dctmtx(size_patch)');
    Height = size(encoded_frames,1);
    Width = size(encoded_frames,2);
    number_of_frames = size(encoded_frames, 3);
    reconstructed_frames = zeros(Height, Width, number_of_frames);
    count_matrix = zeros(Height, Width, number_of_frames);
    for h = 1 : Height - size_patch + 1
        for w = 1: Width - size_patch + 1
            encoding_patch = encoding_code(h:h+size_patch-1, h:h+size_patch-1);
            A = [];
            for i=1:number_of_frames
                ith_frame_patch = encoding_patch(:, :, i);
                phi_t = diag(reshape(ith_frame_patch, size_patch*size_patch, 1));
                a_t = phi_t * dct_2d_psi_matrix;
                A = horzcat(A, a_t);
            end
        frames_patch = encoded_frames(h:h+size_patch-1, h:h+size_patch-1);
        y_vector = reshape(frames_patch, size_patch*size_patch, 1);
        theta = omp_for_loop(A, y, error_limit);

        end
    end
end


    


function theta = omp_for_loop(A, y, error_limit)
    theta = zeros(size(A,2), 1);
    r = y;
    support_set = [];
    Count_index = 0;
    



    





