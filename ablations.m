function [output] = ablations()
%ABLATIONS Summary of this function goes here
%   Detailed explanation goes here
    global W_gap;
    global W_syn;
    global Cap;
    global E_syn;
    global Res;
    
    %% Network Connectivity Matrices of TW Circuit (Adjacency matrix)
    load 'Connectivity_matrix.mat'; % Will Load W_GAP and W_SYN
    W_syn = W_SYN;
    W_gap = W_GAP;
    fid=fopen('resultsAVM.txt','wt');
    stimulus{1} = [1];
    stimulus{2} = [2];
    stimulus{3} = [3];
    fprintf(fid,'FROM ;TO ;HowManyAblated ;OutOf ;FWD_TIME ;BWD_TIME');
    fprintf(fid,'\n');
    [rows, columns] = size(W_syn);
    for i = 1:rows
        for j = 1:columns
            if W_syn(i,j)>0 
                for k = 1:W_syn(i,j)
                    temp = W_syn(i,j);
                    W_syn(i,j) = W_syn(i,j) - k;
                    out = mainTWSimulation(stimulus{1},0);
                    fprintf(fid,'%f %f \n', j);
                    fprintf(fid,';');
                    fprintf(fid,'%f %f \n', i);
                    fprintf(fid,';');
                    fprintf(fid,'%f %f \n', k);
                    fprintf(fid,';');
                    fprintf(fid,'%f %f \n', temp);
                    fprintf(fid,';');                    
                    fprintf(fid,'%f %f \n', out(1));
                    fprintf(fid,';');
                    fprintf(fid,'%f %f \n', out(2));
                    fprintf(fid,';');
                    fprintf(fid,'\n');
%                     disp(out);
                    W_syn(i,j) = temp;
                end
            end
        end
    end
    fprintf(fid,'GAP ;============= ;============= ;============= ;============= ;=============');
    fprintf(fid,'\n');
    [rows, columns] = size(W_gap);
    for i = 1:rows
        for j = 1:i
            if W_gap(i,j)>0 
                for k = 1:W_gap(i,j)
                    temp = W_gap(i,j);
                    W_gap(i,j) = W_gap(i,j) - k;
                    W_gap(j,i) = W_gap(j,i) - k;
                    out = mainTWSimulation(stimulus{1},0);
                    fprintf(fid,'%f %f \n', j);
                    fprintf(fid,';');
                    fprintf(fid,'%f %f \n', i);
                    fprintf(fid,';');
                    fprintf(fid,'%f %f \n', k);
                    fprintf(fid,';');
                    fprintf(fid,'%f %f \n', temp);
                    fprintf(fid,';');                    
                    fprintf(fid,'%f %f \n', out(1));
                    fprintf(fid,';');
                    fprintf(fid,'%f %f \n', out(2));
                    fprintf(fid,';');
                    fprintf(fid,'\n');
%                     disp(out);
                    W_gap(i,j) = temp;
                    W_gap(j,i) = temp;
                end
            end
        end
    end    
fclose(fid);
end

