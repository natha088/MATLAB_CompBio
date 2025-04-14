%% Load Images
a1 = imread('a1.tif');
a2 = imread('a2.tif');
a3 = imread('a3.tif');
a4 = imread('a4.tif');
a5 = imread('a5.tif');
wt1 = imread('wt1.tif');
wt2 = imread('wt2.tif');
wt3 = imread('wt3.tif');
wt4 = imread('wt4.tif');
wt5 = imread('wt5.tif');
wt6 = imread('wt6.tif');

[t,b] = extract_neurite(a1);
[t2,b2] = extract_neurite(a2);
[t3,b3] = extract_neurite(a3);
[t4,b4] = extract_neurite(a4);
[t5,b5] = extract_neurite(a5);
[t6,b6] = extract_neurite(wt1);
[t7,b7] = extract_neurite(wt2);
[t8,b8] = extract_neurite(wt3);
[t9,b9] = extract_neurite(wt4);
[t10,b10] = extract_neurite(wt5);
[t11,b11] = extract_neurite(wt6);
%% Analysis
arcneuro_a = [905.809
576.604
948.335
472.667
846.146];
total_a = [t t2 t3 t4 t5]';
cells_a = [30;29;19;25;21];
r = total_a./arcneuro_a;
ram = mean(r);
ras = std(r);

arcneuro_wt = [258.444
810.778
792.035
388.338
429.332
874.896];
total_wt = [t6 t7 t8 t9 t10 t11]';
cells_wt = [28;36;37;30;28;30];
rwt = total_wt./arcneuro_wt;
rwtm = mean(rwt);
rwts = std(rwt);
[h,p] = ttest2(rwt,r,'Vartype','unequal');
%% Plot Raw & Processed Images
figure
subplot 221
imshow(histeq(imcrop(wt1,[0.5 7.5 1280 820])))
title('1. WT aSyn Raw Image')
subplot 222
imshow(b6)
title('1. WT aSyn Neurite Extraction')
subplot 223
imshow(histeq(imcrop(wt2,[0.5 7.5 1280 820])))
title('2. WT aSyn Raw Image')
subplot 224
imshow(b7)
title('2. WT aSyn Neurite Extraction')

figure
subplot 221
imshow(histeq(imcrop(a1,[0.5 7.5 1280 820])))
title('1. A53T aSyn Raw Image')
subplot 222
imshow(b)
title('1. A53T aSyn Neurite Extraction')
subplot 223
imshow(histeq(imcrop(a2,[0.5 7.5 1280 820])))
title('2. A53T aSyn Raw Image')
subplot 224
imshow(b2)
title('2. A53T aSyn Neurite Extraction')

figure
notBoxPlot([rwt,zeros(6,1)],1:2)
hold on
notBoxPlot([zeros(5,1),r],1:2)
axis([0 3 50 500])
title('Automated Neurosphere Neurite Quantification')
ylabel('Total Dendrite Length per Neurosphere Arc Length')
set(gca,'xtick',[1:2],'xticklabel',{'WT aSyn','A53T aSyn'})
legend(strcat('p-value = ',num2str(p)))


% Sensitivity of mean neurite/neurosphere arc to thresholding:

threshold = [0.20 0.22 0.24];
mean_a53t = [131.7089 167.5836 231.3379];
sd_a53t = [55.7515 49.8584 27.9830];

mean_wt = [243.0426 284.8532 411.1494];
sd_wt = [125.1122 107.869 206.9957];

pvalue = [0.09 0.0478 0.0869];

figure
errorbar(threshold,mean_a53t,sd_a53t)
hold on 
errorbar(threshold,mean_wt,sd_wt)
hold off

figure
imshow(labeloverlay(imcrop(wt2,[0.5 7.5 1280 820]),bwlabel(b7)))

%% ImageJ Analysis Summary

ijrwt = [19643.465
20848.93
18278.85
12777.772
17340.262
14462.182]./arcneuro_wt;

rwtm2 = mean(ijrwt);
rwts2 = std(ijrwt);

ijr = [14994.148
13364.172
12823.397
10857.109
12159.366]./arcneuro_a;

ram2 = mean(ijr);
ras2 = std(ijr);

[h2,p2] = ttest2(ijrwt,ijr,'Vartype','unequal');

figure
notBoxPlot([ijrwt,zeros(6,1)],1:2)
hold on
notBoxPlot([zeros(5,1),ijr],1:2)
axis([0 3 1 100])
title('Manual Neurosphere Neurite Quantification (ImageJ)')
ylabel('Total Dendrite Length per Neurosphere Arc Length')
set(gca,'xtick',[1:2],'xticklabel',{'WT aSyn','A53T aSyn'})
legend(strcat('p-value = ',num2str(p2)))
%% Percent change of A53T to WT 
p_neurite = ((ram - rwtm)/rwtm)*100;
p_ij = ((ram2 - rwtm2)/rwtm2)*100;
error_neurite = p_neurite*sqrt(((ras/ram)^2) + ((rwts/rwtm)^2));
error_ij = p_ij*sqrt(((ras2/ram2)^2) + ((rwts2/rwtm2)^2));

figure
bar([p_neurite,p_ij])
hold on
er = errorbar(1:2,[p_neurite,p_ij],[error_neurite,error_ij],...
    [error_neurite,error_ij]);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 2;
hold off
title({'Percent Change in Neurite Length/Neurosphere Arc Length','between A53T and WT'},'FontSize',14)
xlabel('Quantification Method','FontSize',14)
set(gca,'xtick',[1:2],'xticklabel',{'Automated: Neurite','Manual: ImageJ'},'FontSize',13)
ylabel('% Percent Change','FontSize',14)
%%
function [t,b] = extract_neurite(a1)
a1 = mat2gray(imcrop(sum(a1,3),[0.5 7.5 1280 820]));
a1 = histeq(a1);
a1 = a1 - imopen(a1,strel('disk',20));
a1 = a1 > 0.22;
a1 = imclose(a1,strel('disk',5)) - a1;
a1 = logical(a1) - bwpropfilt(logical(a1), 'Area',[0 50]);
b = zeros(size(a1));
for ii = 0:1:360
    b = b + imopen(a1,strel('line',25,ii));
end
t = sum(sum(a1));
end

