% ECE501B Homework 4
% John Claus
% 10/14/2019

clear all;
close all;


% Define variables
vector_list = [];
final_vector_list = []; 
M = [ 1  0 -1  1  1;
      0  1  0  0  1;
     -1  0 -1 -1 -1;
      1  0 -1  0  1;
      1  1 -1  1 -1;
    ];

%Create a 3D vector list (5x25x5) from 5x25 vector lists
for k = 1:5
    vector_v = rand(5,1); % Radom 5x1 vector
    max_val_v = max(abs(vector_v)); % Used for normaliztion 
    vector_v = vector_v / max_val_v;
%Iterate the input vector into the M matrix 25 times and record on a 5x25 vector list 
    for i = 1:25 
        vector_w = M * vector_v; % Multiply M time the input matrix
        max_val_w = max(abs(vector_w)); % Used for normalization
        vector_w = vector_w / max_val_w; 
        vector_v = vector_w; % Reinject the output into the next interation as an input vector
        vector_list = cat(2, vector_list, vector_w); % Record this iteration's output
    end
    final_vector_list = cat(3, final_vector_list, vector_list); % Add this set of 5x25 vectors to the final list
    vector_list = []; % Clear the old vector list for a new set of 25
end

final_vector_list(:,25,5); % Used to verify the list is correct (troubleshooting)

% Used for part D to find W
[V,D] = eig(M);  
[dummyVar maxEigen] = max(abs(diag(D))); 
W = V(:,maxEigen);  
[dummyVar maxElement] = max(abs(W)); 
W = W ./ W(maxElement); 

% Plots of values calculated earlier
iter = [1 3 5 10 25];
plot(squeeze(final_vector_list(:,iter,5)));
title("Single Input - Various Interation Output");
legend("1 Iteration","3 Iterations","5 Iterations","10 Iterations","25 Iterations");

figure;
iter = [1 2 3 4 5];
plot((squeeze(final_vector_list(:,25,iter)))); 
title("5 Vector Input - 25 Iteration Ouptut");
legend("Vector 1","Vector 2","Vector 3","Vector 4","Vector 5");

figure;
iter = [1 3 5 10 25];
plot(squeeze(final_vector_list(:,iter,5)));
hold on;
plot(W,"r*");
hold off
title("Normalized Eigenvector");
legend("1 Iteration", "3 Iterations", "5 Iterations", "10 Iterations", "25 Iterations","W");

