%% CS754-2022, Assignment 1
% Arpon Basu and Shashwat Garg
% 200050013 & 200050130

% We start with adding the code for the reading the cars
% We extract first T frames, and convert to grayscale


% We obtain the T frames in the first step directly, rather
% than going for T=3, 5 and 7 separately.
file_name='cars.avi';
% file_name='flame.avi';
% Please choose and un-comment the file you wish to use.

addpath('./MMread');
T = 7;

video_data=mmread(file_name,1:T,[],false,true);

H = 120;
W = 240;
epsilon = 48;

Initial_H=video_data.height;
Initial_W=video_data.width;

frame_list=zeros(H,W,T);

for index=1:T   
    frame=rgb2gray(video_data.frames(index).cdata);
    frame_list(:,:,index)=frame(Initial_H-H:Initial_H-1,Initial_W-W:Initial_W-1);
    fig = figure('name',sprintf('Frame %d',index));
    imshow(uint8(frame_list(:,:,index)));
    title(sprintf('Frame %d',index));
    saveas(fig,sprintf('Frame %d.png',index));
end
% By now, we have converted the first 7 frames to grayscale
% We have also extracted the bottom-right section of size 120*240 from all
% the frames


% Generation of Random Code matrix.

Random_C=randi([0 1],H,W,T);


gaussian_noise_stddev = 2;
gaussian_noise_mean = 0;

gaussian_noise=gaussian_noise_mean + rand([H,W])*gaussian_noise_stddev;
% Creation of the gaussian noise matrix to be added to the coded snapshot.


frame_list_3 = frame_list(:,:,1:3);
%This gives us the required working set of frames.
Random_C_3 = Random_C(:,:,1:3);
% This gives us the required measurement matrices which are randomly generated bernoulli matrices. 
E_3 = generate_snapshot(H,W,frame_list_3, 3, Random_C_3,gaussian_noise);
%Coded snapshot has been generated.
I_recon_3 = reconstruction_omp(Random_C_3, E_3, 8, epsilon);
Mean_Squared_Error_3 = mean((I_recon_3 - frame_list_3).^2,'all') / mean(frame_list_3.^2, 'all')




frame_list_5 = frame_list(:,:,1:5);
%This gives us the required working set of frames.
Random_C_5 = Random_C(:,:,1:5);
% This gives us the required measurement matrices which are randomly generated bernoulli matrices. 
E_5 = generate_snapshot(H,W,frame_list_5, 5, Random_C_5,gaussian_noise);
I_recon_5 = reconstruction_omp(Random_C_5, E_5, 8, epsilon);
Mean_Squared_Error_5 = mean((I_recon_5 - frame_list_5).^2,'all') / mean(frame_list_5.^2, 'all')





frame_list_7 = frame_list(:,:,1:7);
%This gives us the required working set of frames.
Random_C_7 = Random_C(:,:,1:7);
% This gives us the required measurement matrices which are randomly generated bernoulli matrices. 
E_7 = generate_snapshot(H,W,frame_list_7, 7, Random_C_7,gaussian_noise);
I_recon_7 = reconstruction_omp(Random_C_7, E_7, 8, epsilon);
Mean_Squared_Error_7 = mean((I_recon_7 - frame_list_7).^2,'all') / mean(frame_list_7.^2, 'all')



















%% Helper functions.

function E = generate_snapshot(Height,Width,frame_list, number_of_frames, Random_C,gaussian_noise)
    E=zeros(Height,Width);
    for i=1:number_of_frames
        E=E+frame_list(:,:,i).*Random_C(:,:,i);
    end
    E = E + gaussian_noise;
    E = max(E, 0);
    E = min(E, number_of_frames*255);
    fig = figure('name',sprintf('Coded snapshot with noise for T=%d',number_of_frames));
    imshow(E/(number_of_frames*255));
    title(sprintf('Coded snapshot with noise for T=%d',number_of_frames));
    saveas(fig,sprintf('Coded_Snapshot_T=%d.png',number_of_frames));

end


function reconstructed_frames = reconstruction_omp(encoding_code, encoded_frames, size_patch, error_limit)
    dct_2d_psi_matrix = kron(dctmtx(size_patch)',dctmtx(size_patch)');
    Height = size(encoding_code,1);
    Width = size(encoding_code,2);
    number_of_frames = size(encoding_code, 3);
    reconstructed_frames = zeros(Height, Width, number_of_frames);
    count_matrix = zeros(Height, Width, number_of_frames);
    square_size = size_patch*size_patch;
    for h = 1 : Height - size_patch + 1
%         h
%         to keep track of algorithm
        for w = 1 : Width - size_patch + 1
            
            encoding_patch = encoding_code(h:h+size_patch-1, w:w+size_patch-1,:);
            A = zeros(square_size, square_size*number_of_frames);
            for t=1:number_of_frames
                ith_frame_patch = encoding_patch(:, :, t);
                phi_t = diag(reshape(ith_frame_patch, square_size, 1));
                a_t = phi_t * dct_2d_psi_matrix;
                A(:,(t-1)*square_size+1:t*square_size) = a_t;
            end
        frames_patch = encoded_frames(h:h+size_patch-1, w:w+size_patch-1);
        y_vector = reshape(frames_patch, square_size, 1);

        theta = omp_for_loop(A, y_vector, error_limit, square_size);
        
        for k=1:number_of_frames
            i_t = dct_2d_psi_matrix * theta((k-1)*size_patch*size_patch+1:k*size_patch*size_patch,:);
            f = reshape(i_t, size_patch, size_patch);
            f = max(f,0);
            f = min(f, 255);
            reconstructed_frames(h:h+size_patch-1, w:w+size_patch-1,k) = reconstructed_frames(h:h+size_patch-1, w:w+size_patch-1,k) + f;
            count_matrix(h:h+size_patch-1, w:w+size_patch-1,k) = count_matrix(h:h+size_patch-1, w:w+size_patch-1,k) + ones(size_patch,size_patch);
        end
        end
    end
    count_matrix = max(count_matrix, 1);
    reconstructed_frames = reconstructed_frames./count_matrix;

    for i=1:number_of_frames
        title_name = sprintf('Reconstruction for T=%d, frame number=%d',number_of_frames, i);
        fig = figure('name',title_name);
        imshow(uint8(reconstructed_frames(:,:,i)));
        title(title_name);
        saveas(fig,sprintf('Reconstruction- T=%d, frame number=%d.png',number_of_frames, i));
    end
end


    


function theta = omp_for_loop(A, y, error_limit, square_size)
    r = y;
    support_set = zeros(1,size(A,2));
    Count_index = 0;
    normalised_A = A./sqrt(sum(A.^2));
    theta = zeros(size(A, 2), 1);
    norm_r = sqrt(sum(r.^2));
    while(norm_r>error_limit && Count_index<square_size)
        Count_index = Count_index+1;
        [~, max_val_index] = max(abs(r' * normalised_A));
        support_set(Count_index) = max_val_index;
        A_support = A(:,support_set(1:Count_index));
        theta(support_set(1:Count_index)) = pinv(A_support)* y;
        r = y - A_support * theta(support_set(1:Count_index));
        norm_r = sqrt(sum(r.^2));
    end
end




      




    



    





