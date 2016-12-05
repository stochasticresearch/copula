% script which generates figure 1& 2 in rsdm paper

fontSize = 20;

M = 500;

figure;
x = rand(M,1);
y = rand(M,1);
u = pobs(x);
v = pobs(y);
p = scatter(u,v);
title('$$X \perp Y$$', 'Interpreter', 'Latex', 'FontSize', fontSize);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
xlabel('$$u = F_X(x)$$', 'FontSize', fontSize, 'Interpreter', 'Latex');
ylabel('$$v = F_Y(y)$$', 'FontSize', fontSize, 'Interpreter', 'Latex');

figure;
x = rand(M,1);
y = sin(2*pi*x);
u = pobs(x);
v = pobs(y);
p = scatter(u,v);
title('$$Y = sin(2 \pi X)$$', 'Interpreter', 'Latex', 'FontSize', fontSize);
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
xlabel('$$u = F_X(x)$$', 'FontSize', fontSize, 'Interpreter', 'Latex');
ylabel('$$v = F_Y(y)$$', 'FontSize', fontSize, 'Interpreter', 'Latex');