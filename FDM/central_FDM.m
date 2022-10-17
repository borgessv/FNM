function [du,O_h] = central_FDM(du_order,n_neighbors)

if mod(n_neighbors,2) ~= 0
    error('The number of neighbors must be an even integer for central finite differences method!')
else
end

h_mult = [n_neighbors/2:-1:1 1:n_neighbors/2];
if mod(du_order,2) == 0
    n_terms = n_neighbors + 1;
else
    n_terms = n_neighbors;
end
syms u [1 n_neighbors+1]
syms du [1 n_terms]
syms alfa [1 n_neighbors]
syms O [1 n_neighbors]
syms h u_j
assume(h,'positive')
u_n = u(n_neighbors/2+1)*sym(ones(n_neighbors,1));
coeff_eq = sym(zeros(n_neighbors,n_terms+1));

for i = 1:n_neighbors
    for j = 1:n_terms
        if i <= n_neighbors/2
            u_n(i) = u_n(i) + 1/factorial(j)*du(j)*(-h*h_mult(i))^(j);
        else
            u_n(i) = u_n(i) + 1/factorial(j)*du(j)*(h*h_mult(i))^(j);
        end
    end
    u_n(i) = u_n(i) + O(i)*h^(j+1);
    coeff_eq(i,:) = coeffs(u_n(i),du);
end

eq_main = sum(u([1:n_neighbors/2 n_neighbors/2+2:end]).*alfa) == sum(u_n.*alfa.');

coeff_eq = sum(coeff_eq(:,2:end),2);
coeff_eq = sum(coeff_eq.*alfa.');
subexp = children(collect(coeff_eq,'h'));
subexp = flip([subexp{:}]);
coeff_sol = solve([subexp(1:du_order-1),subexp(du_order)-1,subexp(du_order+1:end)],alfa);
coeff_sol = struct2cell(coeff_sol);
coeff_sol = [coeff_sol{:}];

sol_main = simplify(solve(subs(eq_main,alfa,coeff_sol),du(du_order)));
O_h = children(collect(sol_main,O),1:n_neighbors);
O_h = coeffs(sum([O_h{:}]),O);
O_h = sum(abs(O_h)./coeffs(sum(abs(O_h)),h));
du = children(collect(sol_main,O),n_neighbors+1);

end
