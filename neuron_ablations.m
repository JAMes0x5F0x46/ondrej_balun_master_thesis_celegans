function [output] = neuron_ablations()
%ABLATIONS Summary of this function goes here
%   Detailed explanation goes here
    global W_gap;
    global W_syn;
    global Cap;
    global E_syn;
    global Res;

    fid=fopen('resultsNeurons.txt','wt');
    stimulus{1} = [1];
    stimulus{2} = [2];
    stimulus{3} = [3];
    stimulus{8} = [8];
    fprintf(fid,'NEURON_ABLATED ;FWD_TIME ;BWD_TIME');
    fprintf(fid,'\n');

    for i = 1:9
        if i ~= 5 && i ~= 7  
            out = mainTWSimulation_neuron_ablation(stimulus{1},i);
            fprintf(fid,'%f %f \n', i);
            fprintf(fid,';');                   
            fprintf(fid,'%f %f \n', out(1));
            fprintf(fid,';');
            fprintf(fid,'%f %f \n', out(2));
            fprintf(fid,';');
            fprintf(fid,'\n');
        end
    end
    
    fprintf(fid,'NEURON_ABLATED ;FWD_TIME ;BWD_TIME');
    fprintf(fid,'\n');

    for i = 1:9
        if i ~= 5 && i ~= 7  
            out = mainTWSimulation_neuron_ablation(stimulus{3},i);
            fprintf(fid,'%f %f \n', i);
            fprintf(fid,';');                   
            fprintf(fid,'%f %f \n', out(1));
            fprintf(fid,';');
            fprintf(fid,'%f %f \n', out(2));
            fprintf(fid,';');
            fprintf(fid,'\n');
        end
    end
fclose(fid);
end

