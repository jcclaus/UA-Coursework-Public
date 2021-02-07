%%%% John Claus %%%%
%%%% ECE501B    %%%%
%%%% Project    %%%%
%%%% 12/08/2019 %%%%

clear all
close all

%American Athletic Conference Matrix 
%Team - Wins - Losses 
AAC = ["Cincinnati",0,0; "UCF",0,0; "Temple",0,0; ...
    "South Florida",0,0; "East Carolina",0,0; "UCONN",0,0; ...
    "Memphis",0,0; "Navy",0,0; "SMU",0,0; ...
    "Tulane",0,0; "Houston",0,0; "Tulsa",0,0];
AAC_length = length(AAC);

%Atlantic Coast Conference Matrix
%Team - Wins - Losses 
ACC = ["Clemson",0,0; "Louisville",0,0; "Wake Forest",0,0; ...
    "Boston College",0,0; "FSU",0,0; "Syracuse",0,0; ...
    "NC State",0,0; "Virginia",0,0; "Virginia Tech",0,0; ...
    "Miami FL",0,0; "North Carolina",0,0; "Pitt",0,0; ...
    "Duke",0,0; "Georgia Tech",0,0];
ACC_length = length(ACC);

%Big 12 Conference Matrix
%Team - Wins - Losses 
B12 = ["Baylor",0,0; "Oklahoma",0,0; "Iowa State",0,0; ...
    "Texas",0,0; "Oklahoma State",0,0; "Kansas State",0,0; ...
    "TCU",0,0; "West Virginia",0,0; "Texas Tech",0,0; ...
    "Kansas",0,0];
B12_length = length(B12);

%Big 10 Conference Matrix
%Team - Wins - Losses 
B10 = ["Wisconsin",0,0; "Minnesota",0,0; "Iowa",0,0; ...
    "Illinois",0,0; "Purdue",0,0; "Nebraska",0,0; ...
    "Northwestern",0,0; "Ohio State",0,0; "Penn State",0,0; ...
    "Michigan",0,0; "Indiana",0,0; "MI State",0,0; ...
    "Maryland",0,0; "Rutgers",0,0];
B10_length = length(B10);

%Conference USA Matrix
%Team - Wins - Losses 
USA = ["UAB",0,0; "LA Tech",0,0; "Southern Miss",0,0; ...
    "Rice",0,0; "North Texas",0,0; "UTSA",0,0; ...
    "UTEP",0,0; "FL Atlantic",0,0; "Marshall",0,0; ...
    "Western Kentucky",0,0; "Charlotte",0,0; "FIU",0,0; ...
    "Middle Tenn",0,0; "Old Dominion",0,0];
USA_length = length(USA);

%Independents FBS Matrix
%Team - Wins - Losses 
IND = ["UMass",0,0; "Army",0,0; "Notre Dame",0,0; ...
    "BYU",0,0; "Liberty",0,0; "NM State",0,0];
IND_length = length(IND);

%Mid-America Conference Matrix
%Team - Wins - Losses 
MAC = ["Miami OH",0,0; "Buffalo",0,0; "Kent State",0,0; ...
    "Bowling Green",0,0; "Akron",0,0; "Central MI",0,0; ...
    "Western MI",0,0; "Ball State",0,0; "Northern IL",0,0; ...
    "Eastern MI",0,0; "Toledo",0,0; "Ohio",0,0];
MAC_length = length(MAC);

%Mountain West Conference Matrix
%Team - Wins - Losses 
MWC = ["Boise State",0,0; "Air Force",0,0; "Utah State",0,0; ...
    "Wyoming",0,0; "CO State",0,0; "New Mexico",0,0; ...
    "Hawaii",0,0; "San Diego St",0,0; "Nevada",0,0; ...
    "UNLV",0,0; "Fresno State",0,0; "San Jose State",0,0];
MWC_length = length(MWC);

%PAC-12 Conference Matrix
%Team - Wins - Losses 
P12 = ["Oregon",0,0; "Cal",0,0; "Washington",0,0; ...
    "Oregon St",0,0; "Stanford",0,0; "Utah",0,0; ...
    "USC",0,0; "UCLA",0,0; "Arizona State",0,0; ...
    "Colorado",0,0; "Arizona",0,0; "Washington State",0,0];
P12_length = length(P12);

%Southeastern Conference Matrix
%Team - Wins - Losses 
SEC = ["LSU",0,0; "Alabama",0,0; "Auburn",0,0; ...
    "Texas AM",0,0; "MS State",0,0; "Miss",0,0; ...
    "Arkansas",0,0; "Georgia",0,0; "Florida",0,0; ...
    "Tennessee",0,0; "South Carolina",0,0; "Kentucky",0,0; ...
    "Missouri",0,0; "Vanderbilt", 0,0];
SEC_length = length(SEC);

%Sun Belt Conference Matrix
%Team - Wins - Losses 
SUN = ["App State",0,0; "GA Southern",0,0; "Georgia State",0,0; ...
    "Troy",0,0; "LA Lafayette",0,0; "Arkansas State",0,0; ...
    "LA Monroe",0,0; "Texas St",0,0; "South Alabama",0,0; ...
    "Coastal",0,0];
SUN_length = length(SUN);


%Import the season records
Records = string(table2cell(readtable('2019_Records.csv')));
rec_length = length(Records);

% Calculate the team Wins and Losses totals
% from the Records data for Colley Method
% games against non-Div1 teams (OTHER) ignored 
for i=1:rec_length
   for j=1:AAC_length
       if (AAC(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               AAC(j,2) = str2num(AAC(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               AAC(j,3) = str2num(AAC(j,3)) + 1;
           end
       end    
   end
   for j=1:ACC_length
       if (ACC(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               ACC(j,2) = str2num(ACC(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               ACC(j,3) = str2num(ACC(j,3)) + 1;
           end
       end    
   end
   for j=1:B12_length
       if (B12(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               B12(j,2) = str2num(B12(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               B12(j,3) = str2num(B12(j,3)) + 1;
           end
       end    
   end
   for j=1:B10_length
       if (B10(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               B10(j,2) = str2num(B10(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               B10(j,3) = str2num(B10(j,3)) + 1;
           end
       end    
   end
   for j=1:P12_length
       if (P12(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               P12(j,2) = str2num(P12(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               P12(j,3) = str2num(P12(j,3)) + 1;
           end
       end    
   end
   for j=1:IND_length
       if (IND(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               IND(j,2) = str2num(IND(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               IND(j,3) = str2num(IND(j,3)) + 1;
           end
       end    
   end
   for j=1:USA_length
       if (USA(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               USA(j,2) = str2num(USA(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               USA(j,3) = str2num(USA(j,3)) + 1;
           end
       end    
   end
   for j=1:MAC_length
       if (MAC(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               MAC(j,2) = str2num(MAC(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               MAC(j,3) = str2num(MAC(j,3)) + 1;
           end
       end    
   end
   for j=1:MWC_length
       if (MWC(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               MWC(j,2) = str2num(MWC(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               MWC(j,3) = str2num(MWC(j,3)) + 1;
           end
       end    
   end
   for j=1:SEC_length
       if (SEC(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               SEC(j,2) = str2num(SEC(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               SEC(j,3) = str2num(SEC(j,3)) + 1;
           end
       end    
   end 
   for j=1:SUN_length
       if (SUN(j,1) == Records(i,1))
           if (Records(i,3) == "W") && (Records(i,2) ~= "OTHER")
               SUN(j,2) = str2num(SUN(j,2)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) ~= "OTHER")
               SUN(j,3) = str2num(SUN(j,3)) + 1;
           end
       end    
   end
end

%Create the Overall Team Matrix
%Team - Wins - Losses 
ALL = [AAC; ACC; B12; B10; USA; IND; MAC; MWC; P12; SEC; SUN];
ALL_length = length(ALL);

% Initialize C and b matrices
C = zeros(length(ALL));
b = zeros(length(ALL),1);

%Create Colley and b matrices
%Create diagnonal values
for j=1:ALL_length
   C(j,j) = str2num(ALL(j,2)) + str2num(ALL(j,3)) + 2;
   b(j) = 1 + ((str2num(ALL(j,2)) - str2num(ALL(j,3)))/2);
end
%Create rest of Colley matrix
for i=1:rec_length
   for j=1:ALL_length
       if (ALL(j,1) == Records(i,1))
          for k=1:ALL_length
            if (ALL(k,1) == Records(i,2))
               C(j,k) = -1; 
            end
          end
       end
   end
end
%Adjustment for Army Navy not yet played 
C(8,66) = 0;
C(66,8) = 0;
%Adjustment for Liberty - NMST game played twice
C(70,69) = -2;
C(69,70) = -2;

%Error check the Colley matrix
%due to large data entry (human error)
%all cells should equal 0 if correct
tot_val = zeros(length(ALL),1);
for i=1:ALL_length
    for j=1:ALL_length
        tot_val(i) = tot_val(i) + C(i,j);
    end
    tot_val(i) = tot_val(i) - 2;
end
error_mat = zeros(length(ALL),2);
index = 1;
%Error checking for matrix structure
%indexes are correclated to team names
for i=1:ALL_length
    for j=1:ALL_length
        if (C(i,j) ~= C(j,i)) && (i < j)
            error_mat(index,1) = i;
            error_mat(index,2) = j;
            index = index + 1;
        end
    end
end
for i=1:ALL_length
   if (error_mat(i,1) > 0) 
       error_teams(i,1) = ALL(error_mat(i,1));
       error_teams(i,2) = ALL(error_mat(i,2));
   end
end

%Solve for the Colley ranking vector
r = linsolve(C,b);

%Create a ranked list of teams for Colley method
Colley_Ranked_Teams = [r+1,ALL(:,1)];
Colley_Ranked_Teams = sortrows(Colley_Ranked_Teams,'descend');
for i=1:length(Colley_Ranked_Teams)
    Colley_Ranked_Teams(i,1) = i;
end

%Create a conference list
conf_list = [0,"AAC"; 0,"ACC"; 0,"B12"; 0,"B10"; 0,"USA"; 0,"IND"; ...
    0,"MAC"; 0,"MWC"; 0,"P12"; 0,"SEC"; 0,"SUN"; 0,"OTH"];

%Import a Colley matrix for Conferences(Claus Matrix)
C_Conf = readmatrix('Claus_Matrix.csv');
b_conf = readmatrix('b_matrix.csv');

%Solve for the conference ranking vector
r_conf = linsolve(C_Conf,b_conf);

%Adjust r vector for non-Div1 conference games
for i=1:11
    r_conf(i) = r_conf(i) - 0.005;
end
r_conf(12) = 11*0.005;

%Combine r vector with conference list
for i=1:length(conf_list)
    conf_list(i,1) = r_conf(i);
end

conf_list = sortrows(conf_list,'descend');

% Calculate the team Wins and Losses totals
% from the Records data for Conference method
% games against non-Div1 teams (OTHER) are not ignored 
CONF_length = ALL_length + 1;
CONF = ALL;
for i=1:CONF_length
    CONF(i,2) = 0;
    CONF(i,3) = 0;
end
CONF(CONF_length,1) = "OTHER";

for i=1:rec_length
   for j=1:ALL_length
       if (CONF(j,1) == Records(i,1))
           if (Records(i,3) == "W") 
               CONF(j,2) = str2num(CONF(j,2)) + 1;
           end
           if (Records(i,3) == "L") 
               CONF(j,3) = str2num(CONF(j,3)) + 1;
           end
           if (Records(i,3) == "W") && (Records(i,2) == "OTHER") 
               CONF(CONF_length,3) = str2num(CONF(CONF_length,3)) + 1;
           end
           if (Records(i,3) == "L") && (Records(i,2) == "OTHER")
               CONF(CONF_length,2) = str2num(CONF(CONF_length,2)) + 1;
           end
       end    
   end
end

% Initialize C and b matrices
C2 = zeros(length(ALL));
b2 = zeros(length(ALL),1);

%Create Colley and b matrices
%Create diagnonal values
for j=1:CONF_length
   C2(j,j) = str2num(CONF(j,2)) + str2num(CONF(j,3)) + 2;
   b2(j) = 1 + ((str2num(CONF(j,2)) - str2num(CONF(j,3)))/2);
end
%Create rest of Colley matrix
for i=1:rec_length
   for j=1:CONF_length
       if (CONF(j,1) == Records(i,1))
          for k=1:CONF_length
            if (CONF(k,1) == Records(i,2))
               C2(j,k) = -1; 
            end
          end
       end
   end
end

%Adjustment for non Div1 games played twice
C2(5,131) = -2;
C2(131,5) = -2;
C2(18,131) = -2;
C2(131,18) = -2;
C2(21,131) = -2;
C2(131,21) = -2;
C2(66,131) = -2;
C2(131,66) = -2;
C2(69,131) = -2;
C2(131,69) = -2;
C2(115,131) = -2;
C2(131,115) = -2;

%Solve for the Colley ranking vector
r2 = linsolve(C2,b2);


for i=1:CONF_length
    if  any(cell2mat(regexp(AAC,CONF(i,1))))
        r2(i) = r2(i) * r_conf(1,1);
    end
    if  any(cell2mat(regexp(ACC,CONF(i,1))))
        r2(i) = r2(i) * r_conf(2,1);
    end
    if  any(cell2mat(regexp(B12,CONF(i,1))))
        r2(i) = r2(i) * r_conf(3,1);
    end
    if  any(cell2mat(regexp(B10,CONF(i,1))))
        r2(i) = r2(i) * r_conf(4,1);
    end
    if  any(cell2mat(regexp(USA,CONF(i,1))))
        r2(i) = r2(i) * r_conf(5,1);
    end
    if  any(cell2mat(regexp(IND,CONF(i,1))))
        r2(i) = r2(i) * r_conf(6,1);
    end
    if  any(cell2mat(regexp(MAC,CONF(i,1))))
        r2(i) = r2(i) * r_conf(7,1);
    end
    if  any(cell2mat(regexp(MWC,CONF(i,1))))
        r2(i) = r2(i) * r_conf(8,1);
    end
    if  any(cell2mat(regexp(P12,CONF(i,1))))
        r2(i) = r2(i) * r_conf(9,1);
    end
    if  any(cell2mat(regexp(SEC,CONF(i,1))))
        r2(i) = r2(i) * r_conf(10,1);
    end
    if  any(cell2mat(regexp(SUN,CONF(i,1))))
        r2(i) = r2(i) * r_conf(10,1);
    end
    if  CONF(i,1) == "OTHER"
        r2(i) = r2(i) * r_conf(11,1);
    end
end


%Create a ranked list of teams for Colley method
Conf_Ranked_Teams = [r2+1,CONF(:,1)];

Conf_Ranked_Teams = sortrows(Conf_Ranked_Teams,'descend');
for i=1:length(Conf_Ranked_Teams)
    Conf_Ranked_Teams(i,1) = i;
end
Ranked_Teams = Colley_Ranked_Teams;
for i=1:length(Colley_Ranked_Teams)
    Ranked_Teams(i,3) = Conf_Ranked_Teams(i,2);
end

