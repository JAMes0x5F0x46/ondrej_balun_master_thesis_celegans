function [output] = VoltageTraces_investor()


 out = simple_investor_Simulation([1], 0);
%out = mainTWSimulation_neuron_ablation(stimulus{1},4)
fid=fopen(strcat('./datasets/investor.txt'),'wt'); 
[r, c] = size(out);
for a = 1:r
    fprintf(fid,'%f %f \n', out(a,1));
    fprintf(fid,';');
    fprintf(fid,'%f %f \n', out(a,2));
    fprintf(fid,';');
    fprintf(fid,'%f %f \n', out(a,3));
    fprintf(fid,';');     
    fprintf(fid,'%f %f \n', out(a,4));
    fprintf(fid,';');   
    fprintf(fid,'%f %f \n', out(a,5));
    fprintf(fid,';');  
    fprintf(fid,'\n');
end

    
end

