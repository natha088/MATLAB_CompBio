% Read DAPI and GFP image files 
u25d = imread('1untreat_20xD.tif'); 
u25g = imread('1untreat_20xG.tif');
u25d2 = imread('1untreat_20x_2D.tif'); 
u25g2 = imread('1untreat_20x_2G.tif');
u25d3 = imread('1untreat_20x_3D.tif'); 
u25g3 = imread('1untreat_20x_3G.tif');
% u25d4 = imread('1untreat_20x_4D.tif'); 
% u25g4 = imread('1untreat_20x_4G.tif');

m25d = imread('1mg132_20xD.tif'); 
m25g = imread('1mg132_20xG.tif');
m25g2 = imread('1mg132_20x_2G.tif');
m25d2 = imread('1mg132_20x_2D.tif');
m25g3 = imread('1mg132_20x_3G.tif');
m25d3 = imread('1mg132_20x_3D.tif');
m25g4 = imread('1mg132_20x_4G.tif');
m25d4 = imread('1mg132_20x_4D.tif');
m25g5 = imread('1mg132_20x_5G.tif');
m25d5 = imread('1mg132_20x_5D.tif');

s25g = imread('1sorb_20xG.tif');
s25d = imread('1sorb_20xD.tif');
s25g2 = imread('1sorb_20x_2G.tif');
s25d2 = imread('1sorb_20x_2D.tif');
s25g3 = imread('1sorb_20x_3G.tif');
s25d3 = imread('1sorb_20x_3D.tif');
s25g4 = imread('1sorb_20x_4G.tif');
s25d4 = imread('1sorb_20x_4D.tif');

% Call CountSG function
[control,d0,b0,g0,dapi0,n0,a0] = CountSG(u25d,u25g);
[control2,d02,b02,g02,dapi02,n02,a02] = CountSG(u25d2,u25g2);
[control3,d03,b03,g03,dapi03,n03,a03] = CountSG(u25d3,u25g3);
% [control4,d04,b04,g04,dapi04,n04,a04] = CountSG(u25d4,u25g4);

[s,d,b,g,dapi,n,a] = CountSG(m25d,m25g);
[s2,d2,b2,g2,dapi2,n2,a2] = CountSG(m25d2,m25g2);
[s3,d3,b3,g3,dapi3,n3,a3] = CountSG(m25d3,m25g3);
[s4,d4,b4,g4,dapi4,n4,a4] = CountSG(m25d4,m25g4);
[s5,d5,b5,g5,dapi5,n5,a5] = CountSG(m25d5,m25g5);

[s10,d10,b10,g10,dapi10,n10,a10] = CountSG(s25d,s25g);
[s12,d12,b12,g12,dapi12,n12,a12] = CountSG(s25d2,s25g2);
[s13,d13,b13,g13,dapi13,n13,a13] = CountSG(s25d3,s25g3);
[s14,d14,b14,g14,dapi14,n14,a14] = CountSG(s25d4,s25g4);

u(:,:,1) = control;
u(:,:,2) = control2;
u(:,:,3) = control3;
% u(:,:,4) = control4;

m(:,:,1) = s;
m(:,:,2) = s2;
m(:,:,3) = s3;
m(:,:,4) = s4;
m(:,:,5) = s5;

s(:,:,1) = s10;
s(:,:,2) = s12;
s(:,:,3) = s13;
s(:,:,4) = s14;
%%
u_rm = mean((u(:,2,:)./u(:,1,:)),3);
m_rm = mean((m(:,2,:)./m(:,1,:)),3);
s_rm = mean((s(:,2,:)./s(:,1,:)),3);

u_rsd = std((u(:,2,:)./u(:,1,:)),[],3);
m_rsd = std((m(:,2,:)./m(:,1,:)),[],3);
s_rsd = std((s(:,2,:)./s(:,1,:)),[],3);

tm = [u_rm,m_rm,s_rm];
ts = [u_rsd,m_rsd,s_rsd];

figure
bar(tm,'grouped')
ax = gca;
ax.FontSize = 16; 
title('Cytoplasmic to Nuclear SG Ratio: 1 ug FUS-GFP','FontSize',18)
xlabel('Threshold (0.132-0.152)','FontSize',18)
ylabel('Cytoplasm:Nuclear SG','FontSize',18)
hold on 
ngroups = size(tm, 1);
nbars = size(tm, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, tm(:,i), ts(:,i), 'k', 'linestyle', 'none');
end
hold off
%% Area

figure
subplot 131
cdfplot([a0.Area]);
hold on;
cdfplot([a02.Area]);
cdfplot([a03.Area]);
hold off
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
title('Untreated: CDF Plot of Stress Granule Area','FontSize',14)
xlabel('Area (pixels)','FontSize',14)
ylabel('Probability','FontSize',14)
subplot 132
cdfplot([a.Area]);
hold on;
cdfplot([a2.Area]);
cdfplot([a3.Area]);
cdfplot([a4.Area]);
cdfplot([a5.Area]);
hold off
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
title('MG-132: CDF Plot of Stress Granule Area','FontSize',14)
xlabel('Area (pixels)','FontSize',14)
ylabel('Probability','FontSize',14)
subplot 133
cdfplot([a10.Area]);
hold on;
cdfplot([a12.Area]);
cdfplot([a13.Area]);
cdfplot([a14.Area]);
hold off
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
title('Sorbitol: CDF Plot of Stress Granule Area','FontSize',14)
xlabel('Area (pixels)','FontSize',14)
ylabel('Probability','FontSize',14)
%% circularity: ratio of minimal to maximal length (r = 1 is circular)
ca0 = [a0.MinorAxisLength]./[a0.MajorAxisLength];
ca02 = [a02.MinorAxisLength]./[a02.MajorAxisLength];
ca03 = [a03.MinorAxisLength]./[a03.MajorAxisLength];

ca = [a.MinorAxisLength]./[a.MajorAxisLength];
ca2 = [a2.MinorAxisLength]./[a2.MajorAxisLength];
ca3 = [a3.MinorAxisLength]./[a3.MajorAxisLength];
ca4 = [a4.MinorAxisLength]./[a4.MajorAxisLength];
ca5 = [a5.MinorAxisLength]./[a5.MajorAxisLength];

ca10 = [a10.MinorAxisLength]./[a10.MajorAxisLength];
ca12 = [a12.MinorAxisLength]./[a12.MajorAxisLength];
ca13 = [a13.MinorAxisLength]./[a13.MajorAxisLength];
ca14 = [a14.MinorAxisLength]./[a14.MajorAxisLength];

circu = [ca0,ca02,ca03];
circm = [ca,ca2,ca3,ca4,ca5];
circs = [ca10,ca12,ca13,ca14];

figure
subplot 131
histfit(circu)
subplot 132 
histfit(circm)
subplot 133
histfit(circs)
%%
% %% Raw Data With Varying Intensity Threshold
% figure
% subplot 231
% bar(control)
% set(gca, 'XTick', 1:11)
% set(gca, 'XTickLabel', 0.2:0.01:0.3)
% xtickangle(60)
% title('Untreated')
% legend('Nuclear','Cytoplasmic','Total')
% xlabel('Intensity Threshold')
% ylabel('Stres Granule Count')
%  
% subplot 232
% bar(s)
% set(gca, 'XTick', 1:11)
% set(gca, 'XTickLabel', 0.2:0.01:0.3)
% xtickangle(60)
% title('MG-132 3')
% legend('Nuclear','Cytoplasmic','Total')
% xlabel('Intensity Threshold')
% ylabel('Stres Granule Count')
%  
% subplot 233
% bar(s1)
% set(gca, 'XTick', 1:11)
% set(gca, 'XTickLabel', 0.2:0.01:0.3)
% xtickangle(60)
% title('Sorbitol')
% legend('Nuclear','Cytoplasmic','Total')
% xlabel('Intensity Threshold')
% ylabel('Stres Granule Count')
%  
% subplot 234
% bar(s2)
% set(gca, 'XTick', 1:11)
% set(gca, 'XTickLabel', 0.2:0.01:0.3)
% xtickangle(60)
% title('Sorbitol 2')
% legend('Nuclear','Cytoplasmic','Total')
% xlabel('Intensity Threshold')
% ylabel('Stres Granule Count')
%  
% subplot 235
% bar(s3)
% set(gca, 'XTick', 1:11)
% set(gca, 'XTickLabel', 0.2:0.01:0.3)
% xtickangle(60)
% title('MG-132')
% legend('Nuclear','Cytoplasmic','Total')
% xlabel('Intensity Threshold')
% ylabel('Stres Granule Count')
%  
% subplot 236
% bar(s4)
% set(gca, 'XTick', 1:11)
% set(gca, 'XTickLabel', 0.2:0.01:0.3)
% xtickangle(60)
% title('MG-132 2')
% legend('Nuclear','Cytoplasmic','Total')
% xlabel('Intensity Threshold')
% ylabel('Stres Granule Count')
% %% Cytoplasmic:Nuclear Stress Granule Ratio with Varying Intensity Threshold
% % Calculate ratio
% ratio(:,1) = control(:,2)./control(:,1);
% ratio(:,2) = s(:,2)./s(:,1);
% ratio(:,3) = s1(:,2)./s1(:,1);
% ratio(:,4) = s2(:,2)./s2(:,1);
% ratio(:,5) = s3(:,2)./s3(:,1);
% ratio(:,6) = s4(:,2)./s4(:,1);
%  
% figure 
% bar(ratio)
% title('Cytoplasmic:Nuclear Stress Granule Ratio','FontSize',14)
% xlabel('Intensity Threshold','FontSize',14)
% ylabel('Cytoplasmic:Nuclear SG Ratio','FontSize',14)
% set(gca, 'XTick', 1:11)
% set(gca, 'XTickLabel', 0.2:0.01:0.3)
% legend('u','m3','s','s2','m','m2')
% %% Area of SG
% 
% figure
% cdfplot([a0.Area]);
% hold on;
% cdfplot([a.Area]);
% cdfplot([a1.Area]);
% cdfplot([a2.Area]);
% cdfplot([a3.Area]);
% cdfplot([a4.Area]);
% hold off
% set(findall(gca, 'Type', 'Line'),'LineWidth',2);
% title('CDF Plot of Stress Granule Area','FontSize',14)
% xlabel('Area (pixels)','FontSize',14)
% ylabel('Probability','FontSize',14)
% legend('u','m3','s','s2','m','m2','FontSize',14)
% 
% %% Pipeline Demo
% figure
% suptitle('Untreated')
% subplot 221
% imshow(imcrop(u25g,[0.5 7.5 1280 820]))
% subplot 222
% imshow(labeloverlay(g0,b0))
% subplot 223
% imshow(imcrop(u25d,[0.5 7.5 1280 820]))
% subplot 224
% imshow(labeloverlay(dapi0,n0))
% 
% figure
% suptitle('Sorbitol')
% subplot 221
% imshow(imcrop(s25g2,[0.5 7.5 1280 820]))
% subplot 222
% imshow(labeloverlay(g2,b2))
% subplot 223
% imshow(imcrop(s25d2,[0.5 7.5 1280 820]))
% subplot 224
% imshow(labeloverlay(dapi2,n2))
% 
% figure
% suptitle('MG-132')
% subplot 221
% imshow(imcrop(m25g,[0.5 7.5 1280 820]))
% subplot 222
% imshow(labeloverlay(g3,b3))
% subplot 223
% imshow(imcrop(m25d,[0.5 7.5 1280 820]))
% subplot 224
% imshow(labeloverlay(dapi3,n3))

%% CountSG function
function [s,count_d,b,g,dapi,nuclei,stats] = CountSG(d,g)
d = mat2gray(imcrop(sum(d,3),[0.5 7.5 1280 820])); % convert images to grayscale and specify dimensions
g = mat2gray(imcrop(sum(g,3),[0.5 7.5 1280 820]));
d = adapthisteq(d);
se = strel('disk',20); % define morphological structure for processing
d = d - imopen(d,se); % opening operation to remove background
g = g - imopen(g,se);
g = imsharpen(g);
d = imbinarize(d);
d = d - bwpropfilt(d,'Area',[0 200]);
d = imdilate(d,strel('disk',1));
dapi = d;
count_d = d;

% Watershed operation
D = bwdist(~d);
D = -D;
L = watershed(D);
count_d(~L) = 0;

count_d = bwlabel(count_d); % label logical 1's
nuclei = count_d; % store total count of nuclei
count_d = max(count_d(:)); % count

s = zeros(11,3);
for ii = 1:11
    g_temp = g > (0.13 + ii/500); % vary intensity threshold to isolate stress granules
    g_temp = bwpropfilt(g_temp,'Area',[5 80]);
    a = bwlabel(d.*g_temp); % find stress granules overlapping with nuclei
    r = bwlabel(g_temp); % find total stress granules
    stats = regionprops('struct',r,'area',...
    'MajorAxisLength','MinorAxisLength');
    if ii == 1
        b = r;
    end
    s(ii,:) = [max(a(:)), max(r(:))-max(a(:)),max(r(:))]; % nuclear, cytoplasmic, total SG
end
end

