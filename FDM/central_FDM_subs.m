function eq = central_FDM_subs(u_var,i,du,h)
eq = subs(du,'h',h);
vars = symvar(eq);
n_var = numel(vars);

if mod(n_var,2) == 0
    neighbors = [-n_var/2:1:-1 1:1:n_var/2];
else
    neighbors = -(n_var-1)/2:1:(n_var-1)/2;
end

eq = subs(eq,vars,u_var(i+neighbors));

end
