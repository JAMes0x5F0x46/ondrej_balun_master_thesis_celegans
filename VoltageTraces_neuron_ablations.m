function [output] = VoltageTraces_neuron_ablations()
%ABLATIONS Summary of this function goes here
%   Detailed explanation goes here
    stimulus{1} = [1];
    stimulus{2} = [2];
    stimulus{3} = [3];
    stimulus{8} = [8];

    for i = 1:9
        if i ~= 5 && i ~= 7  
            out = VoltageTraces_mainTWSimulation_neuron_ablations(stimulus{3},i);
            fid=fopen(strcat('./datasets/neuron_ablation/experiments/',int2str(i),'_plm_missing.txt'),'wt'); 
            [r, c] = size(out);
            for a = 1:r
                fprintf(fid,'%f %f \n', out(a,2));
                fprintf(fid,';');
                fprintf(fid,'%f %f \n', out(a,3));
                fprintf(fid,';');                        
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
    end
end

