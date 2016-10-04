function [result] = nif( LHS, operator, RHS, sigma)

eps = 1e-20;

if  strcmp(operator,'>') 
    diff = LHS - RHS - eps;
    
elseif  strcmp(operator,'>=')
    diff = LHS - RHS;
    
elseif strcmp(operator,'<')
    diff = RHS - LHS - eps;
    
elseif  strcmp(operator,'<=')
    diff = RHS - LHS;
    
elseif strcmp(operator, '==')
    fprintf('Equals are not supported\n')
    return;
    
else
    fprintf('Operator not correct\n')
    return;
end
    p = normcdf(diff,0,sigma);

    q = norminv([(1-p)/2 (1+p)/2], 0, sigma);
    sample = sigma*randn(1,1);

    if (sample >= q(1)) && (sample <= q(2))
        result = true;
    else 
        result = false;
    end

end