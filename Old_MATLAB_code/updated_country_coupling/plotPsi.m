function plotPsi(prop_psi_effects)

prop_psi_effects

fontsize=16;
x_vals = [1 : .1 : 2]
plot(x_vals, prop_psi_effects(:,2), 'r-', x_vals, prop_psi_effects(:,1), 'k-', x_vals, ...
prop_psi_effects(:,3), 'bl-', 'MarkerSize',10, 'LineWidth',2)
legend('Guinea', 'Liberia', 'Sierra Leone', 'Location', 'Northeast')
title('Change in Infected Population Based on Change in \psi','FontSize',fontsize)
xlabel('Multiplicative Increase in \psi','FontSize',fontsize)
ylabel('Infected Proportion','FontSize',fontsize)
set(gca, 'FontSize', fontsize)
print('-painters', '-depsc', '/home/eric/ebola_new/psi_plot.eps')
saveas(gcf,'/home/eric/ebola_new/psi_plot.fig','fig')
