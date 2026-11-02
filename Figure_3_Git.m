close all; clear all ; clc; colors=[[0.97,0.91,0.99];[0.85,0.76,0.89];[0.67,0.53,0.73];[0.51,0.33,0.59];[0.44,0.23,0.53];[0.30,0.09,0.39];[0.5,0.5,0.5]]; set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); fontsize = 20;

rs=[];
Res=logspace(1, 4, 22);
nus=[];

T =4;
nx_pde =500;  nt_pde = 20000; L0 = 1; x0 = -L0/2; x1 = L0/2; t0 = 0; t1 = T; f = 1; U0 = 1/sqrt(nx_pde-1);  

% % %  Grid % % % 
dx_pde = (x1-x0)/(nx_pde-1); dt_pde = (t1-t0)/(nt_pde-1); xs_pde = linspace(x0,x1,nx_pde); ts_pde = linspace(t0,t1,nt_pde);

 % % % Initial conditions % % %

% Original
u0 = @(x) -U0*sin(2*pi*f*x/L0);
% Perturbated
eps = 0.1;
u0_2= @(x) u0(x) - eps*sin(3*2*pi*x/L0);

delta_0= u0(xs_pde)-u0_2(xs_pde);
delta_0_norm= H2norm(xs_pde,delta_0);

% % % Norms % % %
y0_H2=H2norm(xs_pde,u0(xs_pde));


for index=1:length(Res)

    y0_H2= 1.2657;
    Re=Res(index);
    nu=y0_H2/Re
    nus(index)=nu;
    fnu=475.74*nu^2
    fnus=logspace(-6, log10(2*fnu), 16);
    disp(index)

    for index2=1:length(fnus)

        disp(index2)
        fnu_scala=fnus(index2);
        F0_fun = @(t,x) fnu_scala*U0*exp(-(x-L0/4).^2/(2*(L0/32)^2)).*cos(2*pi*t)/ 9.3143;
  
        % % %  PDE, initial condition and boundary condition for pdepe solver % % % 
        burger_pde = @(x,t,u,dudx) deal(1, nu*dudx, -u*dudx + F0_fun(t,x) );
        burger_ic = u0;
        burger_ic_2 =  u0_2;
        burger_bc = @(xl, ul, xr, ur, t) deal(ul, 0, ur, 0);
        
        % % %  Solve exact and perturbed PDEs % % % 
        us_pde = pdepe(0, burger_pde, burger_ic, burger_bc, xs_pde, ts_pde);
        us_pde_2 = pdepe(0, burger_pde, burger_ic_2, burger_bc, xs_pde, ts_pde);
        
        % % % Compute the perturbation over time % % % 
        deltas_norm=[];

        for element=1:length(ts_pde)
            delta=us_pde(element,:)-us_pde_2(element,:);
            deltas_norm(element)=H2norm(xs_pde,delta);
        end
        
        % % % Check if any time the norm of the perturbation was higher than the norm of the initial perturbation % % % 
        rs(index2+16*(index-1))=max(deltas_norm/delta_0_norm);
        
        end
end


%%%%%%%%%%Plotting%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs=fnus/nu^2;
rs_mask=rs;
rs_mask(rs_mask > 1) = 2;
rs_matrix=flipud(reshape(rs_mask, 16, 22));

% Malla de coordenadas
[X, Y] = meshgrid(Res, fs);

% Máscaras
mask1 = (rs_matrix == 1);
mask2 = (rs_matrix == 2);

figure;
hold on;

% Puntos para valores 1
plot(X(mask1), Y(mask1), '.', ...
   'color',colors(5,:), 'MarkerSize', 15);

% Cruces para valores 2
plot(X(mask2), Y(mask2), 'x', ...
   'color',colors(3,:), 'MarkerSize', 8, 'LineWidth', 1.5);

% Ejes personalizados
set(gca, 'XScale', 'log', 'YScale', 'log');


hold on;
maximo_f=1/(3^1.5*0.1/(8*pi^3*0.999984^0.5));
x_curve = logspace(1, log10(max(maximo_f)), 200);
f = @(x) real(4*pi*(0.999984)^0.5/(3^.5*0.100)*cos(pi/3+1/3*acos(3^1.5*0.1*x/(8*pi^3*0.999984^0.5))));  % Define aquí tu función
f2 = @(x) real(4*pi*(0.999984)^0.5/(3^.5*0.100)*cos(pi/3-1/3*acos(3^1.5*0.1*x/(8*pi^3*0.999984^0.5))));  % Define aquí tu función

plot( f(x_curve), x_curve, '-', ...
    'LineWidth', 2, 'color',colors(7,:));

plot( f2(x_curve), x_curve, '-', ...
    'LineWidth', 2, 'color',colors(7,:));

xlim([min(Res), max(1.1*Res)]);
ylim([min(0.9*fs), 1000]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



