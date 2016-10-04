function [output] = VoltageTraces_ablations()
%ABLATIONS Summary of this function goes here
%   Detailed explanation goes here
    global W_gap;
    global W_syn;
    global Cap;
    global E_syn;
    global Res;
    global fid2;
  %  fid2 = fopen(strcat('./datasets/input_current.txt'),'wt');
    
    %% Network Connectivity Matrices of TW Circuit (Adjacency matrix)
    load 'Connectivity_matrix.mat'; % Will Load W_GAP and W_SYN
    W_syn = W_SYN;
    W_gap = W_GAP;
    
    stimulus{1} = [1];
    stimulus{2} = [2];
    stimulus{3} = [3];

    
    %% NO ABLATION DATA
                        out = VoltageTraces_mainTWSimulation(stimulus{1},0);
                    fid=fopen(strcat('./datasets/AVM_stimuli/voltage_traces_no_ablation_AVM_stim_.txt'),'wt'); 
                    [r, c] = size(out);
                    for a = 1:r
                        fprintf(fid,'%f %f \n', out(a,1));
                        fprintf(fid,';');
                        fprintf(fid,'%f %f \n', out(a,2));
                        fprintf(fid,';');
                        fprintf(fid,'%f %f \n', out(a,3));
                        fprintf(fid,';');                        
                        fprintf(fid,'\n');
                    end
                    fclose(fid);
%     
%     [rows, columns] = size(W_syn);
%     for i = 1:rows
%         for j = 1:columns
%             if W_syn(i,j)>0 
%                 for k = 1:W_syn(i,j)
%                     temp = W_syn(i,j);
%                     W_syn(i,j) = W_syn(i,j) - k;
%                     out = VoltageTraces_mainTWSimulation(stimulus{3},0);
%                     fid=fopen(strcat('./datasets/PLM_stimuli/synapses/voltage_traces_s_PLM_stim_',int2str(j),'_',int2str(i),'_',int2str(k),'_',int2str(temp),'.txt'),'wt'); 
%                     [r, c] = size(out);
%                     for a = 1:r
%                         fprintf(fid,'%f %f \n', out(a,1));
%                         fprintf(fid,';');
%                         fprintf(fid,'%f %f \n', out(a,2));
%                         fprintf(fid,';');
%                         fprintf(fid,'%f %f \n', out(a,3));
%                         fprintf(fid,';');                        
%                         fprintf(fid,'\n');
%                     end
%                     W_syn(i,j) = temp;
%                     fclose(fid);
%                 end
%             end
%         end
%     end
%     [rows, columns] = size(W_gap);
%     for i = 1:rows
%         for j = 1:i
%             if W_gap(i,j)>0 
%                 for k = 1:W_gap(i,j)
%                     temp = W_gap(i,j);
%                     W_gap(i,j) = W_gap(i,j) - k;
%                     W_gap(j,i) = W_gap(j,i) - k;
%                     out = VoltageTraces_mainTWSimulation(stimulus{3},0);
%                     fid=fopen(strcat('./datasets/PLM_stimuli/gap_junctions/voltage_traces_g_PLM_stim_',int2str(j),'_',int2str(i),'_',int2str(k),'_',int2str(temp),'.txt'),'wt'); 
%                     [r, c] = size(out);
%                     for a = 1:r
%                         fprintf(fid,'%f %f \n', out(a,1));
%                         fprintf(fid,';');
%                         fprintf(fid,'%f %f \n', out(a,2));
%                         fprintf(fid,';');
%                         fprintf(fid,'%f %f \n', out(a,3));
%                         fprintf(fid,';');                        
%                         fprintf(fid,'\n');
%                     end
%                     W_gap(i,j) = temp;
%                     W_gap(j,i) = temp;
%                     fclose(fid);
%                 end
%             end
%         end
%     end    

% MANUAL ABLATION IN CONNECTIVITY MATRIX

% 
%  out = VoltageTraces_mainTWSimulation(stimulus{3},0);
% %out = mainTWSimulation_neuron_ablation(stimulus{1},4)
% fid=fopen(strcat('./datasets/newTest.txt'),'wt'); 
% [r, c] = size(out);
% for a = 1:r
%     fprintf(fid,'%f %f \n', out(a,1));
%     fprintf(fid,';');
%     fprintf(fid,'%f %f \n', out(a,2));
%     fprintf(fid,';');
% %     fprintf(fid,'%f %f \n', out(a,3));
% %     fprintf(fid,';');     
%     fprintf(fid,'%f %f \n', out(a,4));
%     fprintf(fid,';');
% %     fprintf(fid,'%f %f \n', out(a,5));
% %     fprintf(fid,';');
%     fprintf(fid,'%f %f \n', out(a,6));
%     fprintf(fid,';');  
% %     fprintf(fid,'%f %f \n', out(a,7));
% %     fprintf(fid,';');     
% %     fprintf(fid,'%f %f \n', out(a,8));
% %     fprintf(fid,';');     
%     fprintf(fid,'\n');
% end

    
end

