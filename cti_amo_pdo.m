clear all;

lats = ncread('ersst.v5.185401.nc','lat');
lons = ncread('ersst.v5.185401.nc','lon');
mon{1} = '01'; mon{2} = '02'; mon{3} = '03'; mon{4} = '04'; mon{5} = '05'; mon{6} = '06'; mon{7} = '07'; mon{8} = '08';
mon{9} = '09'; mon{10} = '10'; mon{11} = '11'; mon{12} = '12';
arw = repmat(cosd(double(lats)),[1 numel(lons)]);
arw = arw';
[mask,mask3] = mask_region(lons,lats,sst(:,:,1));

minyear = 1854; maxyear = 2019;
for iy = minyear:maxyear;
    if iy == 2019; io = 8; else; io = 12; end
    for mm = 1:io
    nn = ['ersst.v5.',num2str(iy),mon{mm},'.nc'];
    iyy = (iy - minyear)*12 + mm;
    sst(:,:,iyy) = ncread(nn,'sst');
    end
end

sst(sst < -998) = NaN;

t = month_time([1854:2018]);
t1 = month_time(2019:2019);
t(end+1:end+8) = t1(1:8);

gmst = area_weighted_mean(sst,lons,lats);

%% 1981-2010 mean
tt = find(t >= 1981 & t < 2011);
st1 = sst(:,:,tt);
stt(:,:,:) = nanmean(permute(reshape(st1,[numel(lons),numel(lats),12,numel(st1(1,1,:))/12]),[4,1,2,3]));
stt1 = reshape(repmat(stt,[1,1,1,(numel(sst(1,1,:))-8)/12]),[numel(lons),numel(lats),numel(sst(1,1,:))-8]);
stt1(:,:,end+1:end+8) = stt1(:,:,1:8);

%%
ssta = sst - stt1;
gmsta = area_weighted_mean(ssta,lons,lats);

g1 = permute(repmat(gmsta,[1,numel(lons),numel(lats)]),[2,3,1]);
ssts = ssta - g1;


%% CTI: cti, ctis, reg_cti
fprintf('Computing CTI ...\n');
rlon2 = find(lons >= 180 & lons <= 360-90);
rlat2 = find(lats >= -6 & lats <= 6);

a1 = arw(rlon2,rlat2);
s1 = ssta(rlon2,rlat2,:);
s1 = s1.*repmat(a1,[1 1 size(s1,3)]);
cti = nansum1(reshape(s1,[size(s1,1)*size(s1,2) size(s1,3)]),1)/nansum1(a1(:));
cti = cti(:);


%% AMO (Enfield et al. GRL, 2001): amo3, amo3s, reg_amo3
% ten-year running mean of detrended Atlantic SSTA north of Equator
fprintf('Computing AMO 3 ...\n');
sstar = detrend3d(ssta);
s1 = sstar;
s1(mask3 ~= 3) = nan;
rlats = find(lats < 0);
s1(:,rlats,:) = nan;

amo3 = area_weighted_mean(s1,lons,lats);
amo3 = runmean(amo3,13*12+1,1);


%% PDO: Manuta et al. 1997: lpdo, pdo, reg_pdo
fprintf('Computing PDO ...\n');

s1 = ssts; % Using this after the discussion with Mike.
rlat = find(lats < 20);
s1(mask3(:,:,2:end-1) ~= 2) = nan;
s1(:,rlat,:) = nan;
% weighted by sqrt of cosin of latitude
for j = 1:size(s1,2)
    s1(:,j,:) = s1(:,j,:)*sqrt(cosd(lats(j)));
end
N = 40;
[Lpdo,PCpdo,EOFpdo,stdpcpdo,eigvpdo] = myEOF3dn(s1,N);

rlon20 = find(lons >= 120 & lons <= 130);
rlat20 = find(lats >= 20 & lats <= 40);
d1(:,:) = EOFpdo(:,:,1);
d2 = ones(numel(lons),numel(lats));;
d2(rlon20,rlat20) = 2;
d1(d2 ~= 2) = nan;
d3 = area_weighted_mean(d1,lons,lats);

if d3 < 0
    PCpdo(:,1) = -PCpdo(:,1); EOFpdo(:,:,1) = -EOFpdo(:,:,1);
end
lpdo = Lpdo(1);
pdo = PCpdo(:,1);






