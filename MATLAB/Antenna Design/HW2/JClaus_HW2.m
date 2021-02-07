% John Claus
% ECE584 HW2 

clear all;
format short;

% read from sample_data and put it into a matrix for manipulation
    % dataset1
fid=fopen('dataset1.txt','r');       % read from sample_data.txt to fid
In=fscanf(fid,'%g %g',[2,inf]);      % read formatted file from fid and return the result in matrix form
fclose(fid);                         % close fid
theta_deg1=In(1,:);                  % define 1st column of the matrix as theta_deg
U_data_db1=In(2,:);                  % define 2nd column of the matrix as U_data_db
    % dataset2
fid=fopen('dataset2.txt','r');       % read from sample_data.txt to fid
In=fscanf(fid,'%g %g',[2,inf]);      % read formatted file from fid and return the result in matrix form
fclose(fid);                         % close fid
theta_deg2=In(1,:);                  % define 1st column of the matrix as theta_deg
U_data_db2=In(2,:);                  % define 2nd column of the matrix as U_data_db

% plot the data as two separate figures
    % dataset1
figure(1)                               % define figure #1
plot(theta_deg1,U_data_db1,'b-')        % plot the data in x and y coordinates
xlabel('\theta (degrees)')              % label x coordinate
ylabel('Normalized Power (dB)')         % label y coordinate
title('Normalized Radiation Pattern - dataset 1')   % title of the figure
axis([0 180 -60 0])                     % define the min/max of the x and y axis
    % dataset2
figure(2)                               % define figure #2
plot(theta_deg2,U_data_db2,'b-')        % plot the data in x and y coordinates
xlabel('\theta (degrees)')              % label x coordinate
ylabel('Normalized Power (dB)')         % label y coordinate
title('Normalized Radiation Pattern - dataset 2')   % title of the figure
axis([0 180 -60 0])                     % define the min/max of the x and y axis

% Setup the functions to be numerically integrated
    % dataset1
theta_rad1=theta_deg1*pi/180;                     % convert theta to radians
dtheta1=theta_rad1(2)-theta_rad1(1);              % calculate change in theta per measure 
U_data1=10.^(U_data_db1/10);                      % convert from dB
U_max1=max(U_data1);                              % find max U
F1=U_data1.*sin(theta_rad1)*dtheta1;              % define the integration function
    % dataset2
theta_rad2=theta_deg2*pi/180;                     % convert theta to radians
dtheta2=theta_rad2(2)-theta_rad2(1);              % calculate change in theta per measure 
U_data2=10.^(U_data_db2/10);                      % convert from dB
U_max2=max(U_data2);                              % find max U
F2=U_data2.*sin(theta_rad2)*dtheta2;              % define the integration function

% Calculate Directivity - Piecewise Numerical Integration
    % dataset1
P_rad1=2*pi*sum(F1)         % calculate radiated power
D_01=(4*pi*U_max1)/P_rad1   % calculate directivity
    % dataset2
P_rad2=2*pi*sum(F2)         % calculate radiated power
D_02=(4*pi*U_max2)/P_rad2   % calculate directivity

% Calculate Directivity - Trapezoidal Numerical Integration
    %dataset1
P_rad1T=2*pi*trapz(F1)      % calculate radiated power
D_01T=(4*pi*U_max1)/P_rad1T  % calculate directivity
    % dataset2
P_rad2T=2*pi*trapz(F2)      % calculate radiated power
D_02T=(4*pi*U_max2)/P_rad2T  % calculate directivity


% Creates an Array for Directivity as a function of number of data points
    % dataset1
sz=size(U_data1);                % determines the array size of U_data1
for b = 1:180                    % iterativly calculates the directivity                
    count=0;                     %   by creating new temp arrays composed
    for a = 1:sz(2)              %   of divisor sets of the original
        c=rem(a,b);              % Only records the value into the new temp 
        if c == 0                %    array if the remainder of the divsor 
            count=count+1;       %    value and the current index is zero
            U_data1_Temp(count)=U_data1(a);
            U_data2_Temp(count)=U_data2(a);
            theta_rad1_Temp(count)=theta_rad2(a);
            theta_rad2_Temp(count)=theta_rad2(a);
        end    
    end
    dtheta1_Temp=theta_rad1_Temp(2)-theta_rad1_Temp(1); %calculates the d theta val
    dtheta2_Temp=theta_rad2_Temp(2)-theta_rad2_Temp(1); %calculates the d theta val
    U_max1_Temp=max(U_data1_Temp);                      %finds the new max value of U         
    F1_Temp=U_data1_Temp.*sin(theta_rad1_Temp)*dtheta1_Temp;  %creates the new integration function vals
    U_max2_Temp=max(U_data2_Temp);                      %finds the new max value of U       
    F2_Temp=U_data2_Temp.*sin(theta_rad2_Temp)*dtheta2_Temp;  %creates the new integration function vals
    P1(1,b)=2*pi*sum(F1_Temp);          % First column of power through piecewise
    P1(2,b)=count;                      % Second column of total samples
    P1T(1,b)=2*pi*trapz(F1_Temp);       % First column of power through trapezoidal
    P1T(2,b)=count;                     % Second column of total samples
    P2(1,b)=2*pi*sum(F2_Temp);          % First column of power through piecewise
    P2(2,b)=count;                      % Second column of total samples                     
    P2T(1,b)=2*pi*trapz(F2_Temp);       % First column of power through trapezoidal
    P2T(2,b)=count;                     % Second column of total samples
    D1(1,b)=(4*pi*U_max1_Temp)/P1(1,b); % First column of directivity through piecewise
    D1(2,b)=count;                      % Second column of total samples          
    D1(3,b)=dtheta1_Temp;               % Third column of current delta theta
    D1T(1,b)=(4*pi*U_max1_Temp)/P1T(1,b); % First column of directivity through trapezoidal
    D1T(2,b)=count;                     % Second column of total samples
    D1T(3,b)=dtheta1_Temp;              % Third column of current delta theta
    D2(1,b)=(4*pi*U_max2_Temp)/P2(1,b); % First column of directivity through piecewise  
    D2(2,b)=count;                      % Second column of total samples
    D2(3,b)=dtheta2_Temp;               % Third column of current delta theta
    D2T(1,b)=(4*pi*U_max2_Temp)/P2T(1,b); % First column of directivity through trapezoidal
    D2T(2,b)=count;                     % Second column of total samples
    D2T(3,b)=dtheta2_Temp;              % Third column of current delta theta
    clear *Temp                         % Clears all temp registers
end

% plot the data for directivity as two separate figures
    % dataset1
figure(3)                                             % define figure #3
plot(D1(3,:),D1(1,:),'b-',D1T(3,:),D1T(1,:),'r--')    % plot the data in x and y coordinates
xlabel('d\theta (radians)')                           % label x coordinate
ylabel('Calculated Directivity')                      % label y coordinate
title('Directivity as a Function of d\theta - dataset 1')         % title of the figure
legend('Piecewise','Trapezoidal')                     % differientiates lines on graph
axis([0 1.6 0 7])                                     % define the min/max of the x and y axis
            % samples taken zoomed in at breaking point
figure(4)                                             % define figure #4
plot(D1(2,:),D1(1,:),'b-',D1T(2,:),D1T(1,:),'r--')    % plot the data in x and y coordinates
xlabel('Samples Taken')                               % label x coordinate
ylabel('Calculated Directivity')                      % label y coordinate
title('Directivity as a Function of Samples Taken - dataset 1')   % title of the figure
legend('Piecewise','Trapezoidal')                     % differientiates lines on graph
axis([10 70 0 7])                                     % define the min/max of the x and y axis
    % dataset2
figure(5)                                             % define figure #5
plot(D2(3,:),D2(1,:),'b-',D2T(3,:),D2T(1,:),'r--')    % plot the data in x and y coordinates
xlabel('d\theta (radians)')                           % label x coordinate
ylabel('Calculated Directivity')                      % label y coordinate
title('Directivity as a Function of d\theta - dataset 2')         % title of the figure
legend('Piecewise','Trapezoidal')                     % differientiates lines on graph
axis([0 1.6 0 25])                                    % define the min/max of the x and y axis
            % samples taken zoomed in at breaking point
figure(6)                                             % define figure #6
plot(D2(2,:),D2(1,:),'b-',D2T(2,:),D2T(1,:),'r--')    % plot the data in x and y coordinates
xlabel('Samples Taken')                               % label x coordinate
ylabel('Calculated Directivity')                      % label y coordinate
title('Directivity as a Function of Samples Taken - dataset 2')   % title of the figure
legend('Piecewise','Trapezoidal')                     % differientiates lines on graph
axis([10 70 0 25])                                    % define the min/max of the x and y axis

% Calculates the Side Lobe Angle Values
    % dataset1
pks1=findpeaks(U_data1);        % finds the peak values in the U_data1
pksz=size(pks1);                % determines the size of the array of peak values                        
index_found=0;                  % sets the index found flag to zero
index_count=0;                  % sets the index count to zero
start=1;                        % sets the starting point of search at 1
for i = 1:pksz(2)               % iteratively searches and returns the index values
    for j = start:sz(2)         %   where the peak values exist in U_data1 
        if (index_found==0) && (U_data1(j)==pks1(i)) && (pks1(i)<1)
            index_count=index_count+1;          
            index_val(index_count)=j;  
            index_found=1;
        end
    end
    start=index_val(index_count)+1;
    index_found=0;
end
isz=size(index_val);                    % determines the size of the index_val 
for i = 1:isz(2)                        % correlates the index val with the lobe angle
    Lobes1(1,i)=theta_deg1(index_val(i));
end
clear index*                            % clears the index registers
Lobes1(2,:)=round(10*log10(pks1(pks1<1)));
    % dataset2
pks2=findpeaks(U_data2);        % finds the peak values in the U_data2
pksz=size(pks2);                % determines the size of the array of peak values                        
index_found=0;                  % sets the index found flag to zero
index_count=0;                  % sets the index count to zero
start=1;                        % sets the starting point of search at 1
for i = 1:pksz(2)               % iteratively searches and returns the index values
    for j = start:sz(2)         %   where the peak values exist in U_data1 
        if (index_found==0) && (U_data2(j)==pks2(i)) && (pks2(i)<1)
            index_count=index_count+1;          
            index_val(index_count)=j;  
            index_found=1;
        end
    end
    start=index_val(index_count)+1;
    index_found=0;
end
isz=size(index_val);                    % determines the size of the index_val 
for i = 1:isz(2)                        % correlates the index val with the lobe angle
    Lobes2(1,i)=theta_deg2(index_val(i));
end
clear index*                            % clears the index registers
Lobes2(2,:)=round(10*log10(pks2(pks2<1)));
Lobes1
Lobes2

% Create Graph with Side Lobe Info Overlayed
    % dataset1
figure(7)                               % define figure #7
plot(theta_deg1,U_data_db1,'b-')        % plot the data in x and y coordinates
hold on                                 % allows overlay of scatter plot
scatter(Lobes1(1,:),Lobes1(2,:),'d')    % identifies peaks
d2 = 3;                                 % used to position text
d1 = -3;                                % used to position text
ksz=size(Lobes1);                       % defines size of array for loop
for k = 1:ksz(2)                        % iteratively labels each lobe extrema
    text(Lobes1(1,k)-d2,Lobes1(2,k)+d2,['(' num2str(Lobes1(1,k)) ',' num2str(Lobes1(2,k)) ')'])
    
    %text(Lobes1(1,k)+d1,Lobes1(2,k)+d1,num2str(Lobes1(1,k)))
end
xlabel('\theta (degrees)')              % label x coordinate
ylabel('Normalized Power (dB)')         % label y coordinate
title('Normalized Radiation Pattern - dataset 1')   % title of the figure
axis([0 180 -60 0])                     % define the min/max of the x and y axis
    % dataset2
figure(8)                               % define figure #8
plot(theta_deg2,U_data_db2,'b-')        % plot the data in x and y coordinates
hold on                                 % allows overlay of scatter plot
scatter(Lobes2(1,:),Lobes2(2,:),'d')    % identifies peaks
ksz=size(Lobes2);                       % defines size of array for loop
text(Lobes2(1,1)-10,Lobes2(2,1)+2,['(' num2str(Lobes2(1,1)) ',' num2str(Lobes2(2,1)) ')'])
text(Lobes2(1,2)-10,Lobes2(2,2)-2,['(' num2str(Lobes2(1,2)) ',' num2str(Lobes2(2,2)) ')'])
text(Lobes2(1,3)-12,Lobes2(2,3)+2,['(' num2str(Lobes2(1,3)) ',' num2str(Lobes2(2,3)) ')'])
text(Lobes2(1,4)-8,Lobes2(2,4)-2,['(' num2str(Lobes2(1,4)) ',' num2str(Lobes2(2,4)) ')'])
text(Lobes2(1,5)-10,Lobes2(2,5)+2,['(' num2str(Lobes2(1,5)) ',' num2str(Lobes2(2,5)) ')'])
text(Lobes2(1,6)-2,Lobes2(2,6)-2,['(' num2str(Lobes2(1,6)) ',' num2str(Lobes2(2,6)) ')'])
text(Lobes2(1,7)-10,Lobes2(2,7)+1,['(' num2str(Lobes2(1,7)) ',' num2str(Lobes2(2,7)) ')'])
text(Lobes2(1,8)-8,Lobes2(2,8)+2,['(' num2str(Lobes2(1,8)) ',' num2str(Lobes2(2,8)) ')'])
text(Lobes2(1,9)-8,Lobes2(2,9)+2,['(' num2str(Lobes2(1,9)) ',' num2str(Lobes2(2,9)) ')'])
text(Lobes2(1,10)-8,Lobes2(2,10)+2,['(' num2str(Lobes2(1,10)) ',' num2str(Lobes2(2,10)) ')'])
text(Lobes2(1,11)-6,Lobes2(2,11)+2,['(' num2str(Lobes2(1,11)) ',' num2str(Lobes2(2,11)) ')'])
text(Lobes2(1,12)-8,Lobes2(2,12)+1,['(' num2str(Lobes2(1,12)) ',' num2str(Lobes2(2,12)) ')'])
text(Lobes2(1,13),Lobes2(2,13)+1,['(' num2str(Lobes2(1,13)) ',' num2str(Lobes2(2,13)) ')'])
text(Lobes2(1,14)-12,Lobes2(2,14)-2,['(' num2str(Lobes2(1,14)) ',' num2str(Lobes2(2,14)) ')'])
text(Lobes2(1,15)-6,Lobes2(2,15)-2,['(' num2str(Lobes2(1,15)) ',' num2str(Lobes2(2,15)) ')'])
text(Lobes2(1,16)-6,Lobes2(2,16)+2,['(' num2str(Lobes2(1,16)) ',' num2str(Lobes2(2,16)) ')'])
text(Lobes2(1,17)-8,Lobes2(2,17)-2,['(' num2str(Lobes2(1,17)) ',' num2str(Lobes2(2,17)) ')'])
text(Lobes2(1,18)-8,Lobes2(2,18)+2,['(' num2str(Lobes2(1,18)) ',' num2str(Lobes2(2,18)) ')'])
xlabel('\theta (degrees)')              % label x coordinate
ylabel('Normalized Power (dB)')         % label y coordinate
title('Normalized Radiation Pattern - dataset 2')   % title of the figure
axis([0 180 -60 0])                     % define the min/max of the x and y axis

            
            