figure(1)
fplot('5*exp(-2*x)', [0,2]), hold on
fplot('5*exp(-5*x)', [0,2])
axis([0.6 2 0 5])

figure(2)
fplot('0.09158*exp(2*x)', [0,2]), hold on
fplot('0.000227*exp(5*x)', [0,2])
axis([0 1.4 0 5])

%%
figure(3)
fplot('1510.9852518187017*exp(3149.5*x)', [0,9e-4]), hold on
fplot('1295.5798633789123*exp(4657.5*x)', [0,9e-4])
axis([0 11e-5 1200 1923])

figure(4)
fplot('1923*exp(-1749.4990717548328*x)', [0,0.00096894803548260315]), hold on
fplot('1923*exp(-1763.7125697741442*x)', [0,0.00096894803548260315])
plot([0, 9e-4], [1510, 1510], [0, 9e-4], [1295, 1295])
axis([0 0.00096894803548260315 0 1923])

%%
fplot('2000*exp(-10e8*x)', [0,9e-4]), hold on

%%
xmin = 0;
xmax = 0.00037874682832623622;
Tv = 3315;
Tl = 1923;
Tp1 = 1670.1240925980978;
Tp2 = 1500.1411054766063;
bet1 = 24223.094470442316;
bet2 = 23851.422834040844;

X = linspace(xmin,xmax,100);
Tx1 = (Tv-Tp1)*exp(-bet1.*X) + Tp1;
Tx2 = (Tv-Tp2)*exp(-bet2.*X) + Tp2;

figure(3)
plot(X, Tx1, X, Tx2), hold on

%%
xmin = 0;
xmax = 0.0100;
% xmax = 0.00037874682832623622;

Tv = 3315;
Tl = 1923;
Tp1 = 353.1240925980978;
Tp2 = 353.01411054766063;

%----- Original Model -----
bet1 = expCoeff_DELETEME(xmax - xmin, Tp1, Tl);
bet2 = expCoeff_DELETEME(xmax - xmin, Tp2, Tl);

X = linspace(xmin,xmax,100);
Tx1 = Tl*exp(-bet1.*X);
Tx2 = Tl*exp(-bet2.*X);

%----- For Cooling -----
Tmax1 = Tl;
Tmax2 = 1850;
bet1_cool = expCoeff_DELETEME(xmax - xmin, Tp1, Tmax1);
bet2_cool = expCoeff_DELETEME(xmax - xmin, Tp2, Tmax2);

Tx1_cool = Tmax1*exp(-bet1_cool.*X);
Tx2_cool = Tmax2*exp(-bet2_cool.*X);



%----- Plot Resulting Curves ----- 
figure(3)
plot(X, Tx1, X, Tx2, X, Tx1_cool, X, Tx2_cool, '-k'), hold on
legend('Tx1', 'Tx2', 'Tx1_c_o_o_l', 'Tx2_c_o_o_l')

