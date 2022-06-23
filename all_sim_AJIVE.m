
clear
clc

% set path having AJIVE matlab code first
addpath ../AJIVEcode/


ioutput = [1, 1, 1, 1, 1, 1, 1, 1, 1];
paramstruct0 = struct('ioutput', ioutput, ...
                      'iplot', [0 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 1. D = 2 with orthogonal case in the main paper %%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = './simulations/main/2views/orthogonal/dat/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
error = zeros(100,1);


for i=1:100
    
    % Load the data
    filename = [folder 'dat' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = Z1';
    datablock{2} = Z2'; 
    p1 = size(Z1,2);
    p2 = size(Z2,2);

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);

    % Get back the estimation error
    
    for d=1:2
        estimate = fb_vec(d)*(outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d})';
        if d==1
            index = 1:p1;
        else
            index = (p1+1):(p1+p2);
        end
        error(i) = error(i) + sum(sum((true_M(:,index) - estimate).^2))/sum(sum(true_M(:,index).^2));
    end    
    
end

savefile = [folder 'ajive_res.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 2. D = 2 with non-orthogonal case in the main paper %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = './simulations/main/2views/non_orthogonal/dat/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
error = zeros(100,1);


for i=1:100
    
    % Load the data
    filename = [folder 'dat' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = Z1';
    datablock{2} = Z2'; 
    p1 = size(Z1,2);
    p2 = size(Z2,2);

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);

    % Get back the estimation error
    
    for d=1:2
        estimate = fb_vec(d)*(outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d})';
        if d==1
            index = 1:p1;
        else
            index = (p1+1):(p1+p2);
        end
        error(i) = error(i) + sum(sum((true_M(:,index) - estimate).^2))/sum(sum(true_M(:,index).^2));
    end    
    
end

savefile = [folder 'ajive_res.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 3. D = 3 with orthogonal case in the main paper %%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = './simulations/main/3views/orthogonal/dat/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
r3vec = zeros(100,1);
error = zeros(100,1);


for i=1:100
    
    % Load the data
    filename = [folder 'dat' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = Z1';
    datablock{2} = Z2'; 
    datablock{3} = Z3';     
    p1 = size(Z1,2);
    p2 = size(Z2,2);
    p3 = size(Z3,2);

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);
    r3vec(i) = size(outstruct.BSSindiv{3},1);

    % Get back the estimation error
    
    for d=1:3
        estimate = fb_vec(d)*(outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d})';
        if d==1
            index = 1:p1;
        elseif d==2
            index = (p1+1):(p1+p2);
        else
            index = (p1+p2+1):(p1+p2+p3);
        end
        error(i) = error(i) + sum(sum((true_M(:,index) - estimate).^2))/sum(sum(true_M(:,index).^2));
    end    
    
end

savefile = [folder 'ajive_res.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 4. D = 3 with non-orthogonal case in the main paper %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = './simulations/main/3views/non_orthogonal/dat/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
r3vec = zeros(100,1);
error = zeros(100,1);


for i=1:100
    
    % Load the data
    filename = [folder 'dat' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = Z1';
    datablock{2} = Z2'; 
    datablock{3} = Z3';     
    p1 = size(Z1,2);
    p2 = size(Z2,2);
    p3 = size(Z3,2);

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);
    r3vec(i) = size(outstruct.BSSindiv{3},1);

    % Get back the estimation error
    
    for d=1:3
        estimate = fb_vec(d)*(outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d})';
        if d==1
            index = 1:p1;
        elseif d==2
            index = (p1+1):(p1+p2);
        else
            index = (p1+p2+1):(p1+p2+p3);
        end
        error(i) = error(i) + sum(sum((true_M(:,index) - estimate).^2))/sum(sum(true_M(:,index).^2));
    end    
    
end

savefile = [folder 'ajive_res.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 5. 1st rank scheme with orthogonal case in Appendix %%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = './simulations/appendix/setup1/orthogonal/dat/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
r3vec = zeros(100,1);
error = zeros(100,1);


for i=1:100
    
    % Load the data
    filename = [folder 'dat' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = Z1';
    datablock{2} = Z2'; 
    datablock{3} = Z3';     
    p1 = size(Z1,2);
    p2 = size(Z2,2);
    p3 = size(Z3,2);

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);
    r3vec(i) = size(outstruct.BSSindiv{3},1);

    % Get back the estimation error
    
    for d=1:3
        estimate = fb_vec(d)*(outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d})';
        if d==1
            index = 1:p1;
        elseif d==2
            index = (p1+1):(p1+p2);
        else
            index = (p1+p2+1):(p1+p2+p3);
        end
        error(i) = error(i) + sum(sum((true_M(:,index) - estimate).^2))/sum(sum(true_M(:,index).^2));
    end    
    
end

savefile = [folder 'ajive_res.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 5. 1st rank scheme with non-orthogonal case in Appendix %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = './simulations/appendix/setup1/non_orthogonal/dat/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
r3vec = zeros(100,1);
error = zeros(100,1);


for i=1:100
    
    % Load the data
    filename = [folder 'dat' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = Z1';
    datablock{2} = Z2'; 
    datablock{3} = Z3';     
    p1 = size(Z1,2);
    p2 = size(Z2,2);
    p3 = size(Z3,2);

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);
    r3vec(i) = size(outstruct.BSSindiv{3},1);

    % Get back the estimation error
    
    for d=1:3
        estimate = fb_vec(d)*(outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d})';
        if d==1
            index = 1:p1;
        elseif d==2
            index = (p1+1):(p1+p2);
        else
            index = (p1+p2+1):(p1+p2+p3);
        end
        error(i) = error(i) + sum(sum((true_M(:,index) - estimate).^2))/sum(sum(true_M(:,index).^2));
    end    
    
end

savefile = [folder 'ajive_res.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 7. 2nd rank scheme with orthogonal case in Appendix %%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = './simulations/appendix/setup2/orthogonal/dat/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
r3vec = zeros(100,1);
error = zeros(100,1);


for i=1:100
    
    % Load the data
    filename = [folder 'dat' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = Z1';
    datablock{2} = Z2'; 
    datablock{3} = Z3';     
    p1 = size(Z1,2);
    p2 = size(Z2,2);
    p3 = size(Z3,2);

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);
    r3vec(i) = size(outstruct.BSSindiv{3},1);

    % Get back the estimation error
    
    for d=1:3
        estimate = fb_vec(d)*(outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d})';
        if d==1
            index = 1:p1;
        elseif d==2
            index = (p1+1):(p1+p2);
        else
            index = (p1+p2+1):(p1+p2+p3);
        end
        error(i) = error(i) + sum(sum((true_M(:,index) - estimate).^2))/sum(sum(true_M(:,index).^2));
    end    
    
end

savefile = [folder 'ajive_res.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% 8. 2nd rank scheme with non-orthogonal case in Appendix %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = './simulations/appendix/setup2/non_orthogonal/dat/';

r0vec = zeros(100,1);
r1vec = zeros(100,1);
r2vec = zeros(100,1);
r3vec = zeros(100,1);
error = zeros(100,1);


for i=1:100
    
    % Load the data
    filename = [folder 'dat' int2str(i) '.mat'];
    load (filename);

    % Apply AJIVE
    datablock{1} = Z1';
    datablock{2} = Z2'; 
    datablock{3} = Z3';     
    p1 = size(Z1,2);
    p2 = size(Z2,2);
    p3 = size(Z3,2);

    outstruct = AJIVEMainMJ(datablock, vecr, paramstruct0);
    
    % Get back the ranks
    r0vec(i) = size(outstruct.CNS,1);
    r1vec(i) = size(outstruct.BSSindiv{1},1);
    r2vec(i) = size(outstruct.BSSindiv{2},1);
    r3vec(i) = size(outstruct.BSSindiv{3},1);

    % Get back the estimation error
    
    for d=1:3
        estimate = fb_vec(d)*(outstruct.MatrixJoint{d}+ outstruct.MatrixIndiv{d})';
        if d==1
            index = 1:p1;
        elseif d==2
            index = (p1+1):(p1+p2);
        else
            index = (p1+p2+1):(p1+p2+p3);
        end
        error(i) = error(i) + sum(sum((true_M(:,index) - estimate).^2))/sum(sum(true_M(:,index).^2));
    end    
    
end

savefile = [folder 'ajive_res.mat'];
save(savefile, 'r0vec', 'r1vec', 'r2vec', 'error');





