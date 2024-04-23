% Function inputs
SuppPlot = true;
inputCoords = OUTMAT_permuted; %name of your input data?
nmodes = 3; % How many modes do you want to display?
oxDNA_Sim = false;
mapped_crd = false;
sz = size(inputCoords);
npart = sz(1);
range_to_pca = 1:sz(2);
frames_to_remove = [];
remove_last_frame = false;
coords = flip(inputCoords(:,:,:));
coords(:,frames_to_remove,:)=[];
fprintf("Computing mean structure\n");
MeanStruct = computeMeanStructure(coords);
fprintf("Finished computing mean structure\n");
fprintf("Computing deviations from mean structure via procrustes\n");
sz = size(coords);
lsim = sz(2);
flucts = zeros(lsim, npart, 3);
for i=1:lsim
    [d,Z] = procrustes(MeanStruct,squeeze(coords(:,i,1:3)), 'scaling', false, 'reflection', false);
    flucts(i,:,:) = Z-MeanStruct;
    if (mod(i,100)==0)
        fprintf("Fluctuation Calculations %i/%i\n",i,lsim);
    end
end

fprintf("Finished computing structure fluctuations\n");

% Now compute covariance matrix and diagonalize
%fprintf("Computing mean deviation\n");
%MeanDeviation = squeeze(mean(flucts));

Flucts_vec = reshape(flucts,lsim, npart*3);
C = zeros(npart*3, npart*3);     
for i=1:lsim
    C = C + transpose(Flucts_vec(i,:))*Flucts_vec(i,:);
end
C = C / lsim;
e = eig(C);
totalF = sum(e);
FracF = e./totalF;
PctF = FracF * 100;
figure(1)
hold on
scatter(1:length(PctF),flip(PctF), 256, 'o');
xlabel("Principal Component");
ylabel("\lambda/\Sigma_i\lambda_i * 100%")
title("PCA Scree plots, oxDNA and coarsened representation");
plt=Plot();
[V,D] = eig(C);
plt.BoxDim=[12 6];
xlim([0 50]);

figure(4)
% setup color array
carray=[];
for i=1:npart
    carray = [carray; [1,0,0]];
end

for i=1:nmodes
    if (SuppPlot)
        subplot(1,nmodes,i);
    else
        figure(i+1);
    end
    hold off
    XVALS = MeanStruct(:,1);
    YVALS = MeanStruct(:,2);
    ZVALS = MeanStruct(:,3);
    PCAVec = V(:,end-i+1);
    PCAVec = reshape(PCAVec,npart,3);
    % Show the scaffold to have arrows plotted over it
    %plot3(XVALS,YVALS,ZVALS,'LineWidth',12);
    tubeplot(transpose(MeanStruct),.8,16, [], [0.3 0.3 0.8]);
    PCAX = PCAVec(:,1)*e(end-i+1);
    PCAY = PCAVec(:,2)*e(end-i+1);
    PCAZ = PCAVec(:,3)*e(end-i+1);
    hold on
    ScaleDown = 6;
    %NewTubePlotCoords = transpose(MeanStruct)+transpose([PCAX./ScaleDown, PCAY./ScaleDown, PCAZ./ScaleDown]);
    %tubeplot(NewTubePlotCoords,1,16, [], [0 1 1]);
    %NewTubePlotCoords = transpose(MeanStruct)+transpose([-PCAX./ScaleDown, -PCAY./ScaleDown, -PCAZ./ScaleDown]);
    %tubeplot(NewTubePlotCoords,1,16, [], [0 0.5 0.25]);
    quiver3D([XVALS,YVALS,ZVALS],-[PCAX,PCAY,PCAZ]./ScaleDown, carray);
    if (~SuppPlot)
        titlestring = sprintf("Principal Component %i", i);
        title(titlestring);
    end
    
    grid off
    axis off
    set (gcf, 'color','w');
    camlight
    lightangle(-45,30)
end

if (SuppPlot)
    grouptitlestr = sprintf("First %i components of motion", nmodes);
    sgtitle(grouptitlestr);
end