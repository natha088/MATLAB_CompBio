sample_data = xlsread('autocorr_sample_data.xlsx',1,'A1:A30');
t = linspace(0,116,30)';
[a,b,c] = fit_int_autocorr(t,sample_data);
lifetime = c;
rsq = b.rsquare;
figure
scatter(t,sample_data)
hold on;
plot(a)
txt = ['Hydrogen Bond Lifetime: ' num2str(lifetime) ' ps'];
text(50,0.5,txt)
hold off;
xlabel(' t (ps)')
ylabel('autocorr(t)')


