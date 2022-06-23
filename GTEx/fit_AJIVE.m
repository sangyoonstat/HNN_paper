
clear
clc
addpath ../AJIVECode/

load("./GTEx/GTEx_data.mat")

Z1 = muscle/norm(muscle, "fro");
Z2 = blood/norm(blood, "fro");
Z3 = skin/norm(skin, "fro");


ioutput = [1, 1, 1, 1, 1, 1, 1, 1, 1];
paramstruct0 = struct('ioutput', ioutput, ...
                      'iplot', [0 0]);


% Apply AJIVE
 datablock{1} = Z1';
 datablock{2} = Z2'; 
 datablock{3} = Z3';     
 p1 = size(Z1,2);
 p2 = size(Z2,2);
 p3 = size(Z3,2);    
 vecr=[12,5,11];

outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);

% Apply AJIVE to 1st and 2nd views
datablock12{1} = Z1';
datablock12{2} = Z2'; 
vecr12=[12,5];
outstruct12 = AJIVEMainMJ(datablock12, vecr12, paramstruct0);    

% Apply AJIVE to 1st and 3rd views
datablock13{1} = Z1';
datablock13{2} = Z3'; 
vecr13=[12,11];
outstruct13 = AJIVEMainMJ(datablock13, vecr13, paramstruct0);

% Apply AJIVE to 2nd and 3rd views
datablock23{1} = Z2';
datablock23{2} = Z3'; 
vecr23=[5,11];
outstruct23 = AJIVEMainMJ(datablock23, vecr23, paramstruct0);

