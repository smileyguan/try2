%% Generate New Response List by using maximum value between 13-16 and 14-17
%% Usining the temporal data with 24 frames
%% Edited by Guan 20170921

clear;clc

load temporal_tuning.mat;


%% I am a fool %%

%% lala %%

FrameNum = size(temporal_tuning,1);
SfNum = 6;
OriNum = 12;
SzNum = 3;
CellNum = size(temporal_tuning,5);
RepeatNum = size(temporal_tuning,6);

G4_RspAvgOFFListTotal = zeros(SfNum*SzNum,CellNum,OriNum,RepeatNum);


for ci=1:CellNum
    for a = 1:OriNum
        for b = 1:SfNum
            for c = 1:SzNum
                OffImage = mean(mean(squeeze(temporal_tuning(5:8,b,a,c,ci,:)))); %%%% Average 10 trials 5:8 frames as OFFimage for each condition
                for d = 1:RepeatNum
                    OnImage_two = [mean(squeeze(temporal_tuning(13:16,b,a,c,ci,d))),mean(squeeze(temporal_tuning(14:17,b,a,c,ci,d)))];
                    OnImage = max(OnImage_two); %%%% Pick the maximum response between 13:16 frames and 14:17 frames 
                    rsp(d) = (OnImage-OffImage)/OffImage;
                end
                G4_RspAvgOFFListTotal(3*(b-1)+c,ci,a,:) = rsp;
            end
        end
    end
end

save G4_RspAvgOFFListTotal.mat G4_RspAvgOFFListTotal;


%% Pick the maximum sf and orientation
%% Save the mean response of each condition by avg repeat trials
load G4_RspAvgOFFListTotal.mat;



SfSzNum = size(G4_RspAvgOFFListTotal,1);
cellNum = size(G4_RspAvgOFFListTotal,2);
oln = size(G4_RspAvgOFFListTotal,3);
ol = -45:15:120;
trialNum = size(G4_RspAvgOFFListTotal,4);
G4_RspMeanTrialStdSeListTotal = zeros(SfSzNum,cellNum,oln,4);
%1.Useful trial avg(<3sd); 2.Useful trial number; 3.Useful trial std; 4.Std/sqrt(trial number);

for ci = 1:cellNum  %%CellNum
    for ss = 1:SfSzNum
        for j3 = 1:oln %%Orientaion
            rspTemporal = squeeze(G4_RspAvgOFFListTotal(ss,ci,j3,:));
            rspMean = mean(rspTemporal);
            rspStd = std(rspTemporal);
            index = find(rspTemporal>rspMean+3*rspStd);
            index = [index,find(rspTemporal<rspMean-3*rspStd)];
%             rspTemporal(index) = [];     %%%%Dlete >3std trials
            rspTrial = length(rspTemporal);
            rspMean = mean(rspTemporal);
            rspStd = std(rspTemporal);
            rspSe = rspStd./sqrt(rspTrial);
            G4_RspMeanTrialStdSeListTotal(ss,ci,j3,1) = rspMean;
            G4_RspMeanTrialStdSeListTotal(ss,ci,j3,2) = rspTrial;
            G4_RspMeanTrialStdSeListTotal(ss,ci,j3,3) = rspStd;
            G4_RspMeanTrialStdSeListTotal(ss,ci,j3,4) = rspSe;
        end
    end
end

save G4_RspMeanTrialStdSeListTotal.mat G4_RspMeanTrialStdSeListTotal; %%Normal distribution from -45:15:120


%% Find the Maximum Response in the matrix
%%%Pick the maximum response Size in each Sf*Ori condition(3to1)

G4_SfOri_MeanSeListTotal = zeros(SfSzNum/3,cellNum,oln,2);
%1.Mean; 2.Se
G4_SzLocListTotal = zeros(SfSzNum/3,cellNum,oln);


for ci = 1:cellNum
    SfSzNum_Oln_Mean = squeeze(G4_RspMeanTrialStdSeListTotal(:,ci,:,1));
    SfSzNum_Oln_Se = squeeze(G4_RspMeanTrialStdSeListTotal(:,ci,:,4));
    for j3 = 1:oln
        for j4 = 1:SfSzNum/3
            [MaxSzRsp,MaxSz] = max(squeeze(SfSzNum_Oln_Mean(j4*3-2:j4*3,j3)));
            MaxSzSe = SfSzNum_Oln_Se(j4*3-3+MaxSz,j3);
            G4_SfOri_MeanSeListTotal(j4,ci,j3,1) = MaxSzRsp;
            G4_SfOri_MeanSeListTotal(j4,ci,j3,2) = MaxSzSe;
            G4_SzLocListTotal(j4,ci,j3) = j4*3-3+MaxSz;
        end
    end
end


save G4_SfOri_MeanSeListTotal.mat G4_SfOri_MeanSeListTotal;
save G4_SzLocListTotal.mat G4_SzLocListTotal;

%%%Pick the maximum Ori&Sf response

load G4_SfOri_MeanSeListTotal.mat;
load G4_SzLocListTotal.mat

G4_RawDataforFITTINGListTotal = zeros(cellNum,18,2);
G4_PeakSfLocListTotal = zeros(cellNum,20);
%%The 1:12 is different Ori with Max Sf, 13:18 is different Sf with Max Ori


for ci = 1:cellNum
    SfOri_Mean = squeeze(G4_SfOri_MeanSeListTotal(:,ci,:,1));
    SfOri_Se = squeeze(G4_SfOri_MeanSeListTotal(:,ci,:,2));
    SzLoc = squeeze(G4_SzLocListTotal(:,ci,:));
    [PeakSf PeakOri]=find(SfOri_Mean==max(max(SfOri_Mean)));
    G4_RawDataforFITTINGListTotal(ci,1:12,1) = SfOri_Mean(PeakSf,:);
    G4_RawDataforFITTINGListTotal(ci,13:18,1) = (SfOri_Mean(:,PeakOri))';
    G4_RawDataforFITTINGListTotal(ci,1:12,2) = SfOri_Se(PeakSf,:);
    G4_RawDataforFITTINGListTotal(ci,13:18,2) = (SfOri_Se(:,PeakOri))';
    G4_PeakSfLocListTotal(ci,1:12) = SzLoc(PeakSf,:);
    G4_PeakSfLocListTotal(ci,13) = PeakOri; %%13 is peak orientation
    G4_PeakSfLocListTotal(ci,14:19) = SzLoc(:,PeakOri);
    G4_PeakSfLocListTotal(ci,20) = PeakSf;
end


save G4_RawDataforFITTINGListTotal.mat G4_RawDataforFITTINGListTotal;
save G4_PeakSfLocListTotal.mat G4_PeakSfLocListTotal;


%% Find Anova Data

load G4_RspAvgOFFListTotal.mat;
load G4_PeakSfLocListTotal.mat;

cellNum = size(G4_RspAvgOFFListTotal,2);
oln = size(G4_RspAvgOFFListTotal,3);
sfn = 6;
ol = -45:15:120;
trialNum = size(G4_RspAvgOFFListTotal,4);
Y1_AnovaListTotal = zeros(1,cellNum);

for ci = 1:cellNum
    anovaList = zeros(oln+sfn,trialNum);
    for j3 = 1:oln
        Loc = G4_PeakSfLocListTotal(ci,j3);
        anovaList(j3,:)= G4_RspAvgOFFListTotal(Loc,ci,j3,:);
    end
    for j = oln+1:oln+sfn
        Loc2 = G4_PeakSfLocListTotal(ci,j+1);
        PeakOri = G4_PeakSfLocListTotal(ci,13);
        anovaList(j,:)= G4_RspAvgOFFListTotal(Loc2,ci,PeakOri,:);
    end    
        
    anovaList = anovaList';
    anova_test = kruskalwallis(anovaList(:,1:12),1:12,'off');
    p = anova_test(1);
    Y1_AnovaListTotal(1,ci) = p; %%%OriAnova
    
    anova_test2 = kruskalwallis(anovaList(:,13:18),1:6,'off');
    p2 = anova_test2(1);
    F4_AnovaListTotal(1,ci) = p2; %%%SFAnova
end

save Y1_AnovaListTotal.mat Y1_AnovaListTotal;
save F4_AnovaListTotal.mat F4_AnovaListTotal;

%% direction pick
clear;clc

load G4_PeakSfLocListTotal.mat;
load G4_RspAvgOFFListTotal.mat;

Dir_Matrix = zeros(size(G4_PeakSfLocListTotal,1), 12);

for i = 1:size(G4_PeakSfLocListTotal,1)
    ci = i;
    PeakOri = G4_PeakSfLocListTotal(ci, 13);
    PeakSf = G4_PeakSfLocListTotal(ci, 20);
    PeakNum = G4_PeakSfLocListTotal(ci, 13+PeakSf);
    Dir_Matrix(ci, :) = squeeze(G4_RspAvgOFFListTotal(PeakNum, ci, PeakOri, :));
end

save Dir_Matrix Dir_Matrix;

pvalueMat = zeros(size(Dir_Matrix,1),1);
for i = 1:size(Dir_Matrix,1)
    p = ranksum(Dir_Matrix(i,1:6), Dir_Matrix(i,7:12));
    pvalueMat(i) = p;
end

save pvalueMat pvalueMat;


load G4_RspMeanTrialStdSeListTotal;
withdir_index = find(pvalueMat<0.05);

    
SfSzNum = size(G4_RspAvgOFFListTotal,1);
cellNum = size(G4_RspAvgOFFListTotal,2);
oln = size(G4_RspAvgOFFListTotal,3);
ol = -45:15:120;
trialNum = 6;
%1.Useful trial avg(<3sd); 2.Useful trial number; 3.Useful trial std; 4.Std/sqrt(trial number);
for k = 1:length(withdir_index)
    ci = withdir_index(k);
    dirM = [mean(Dir_Matrix(ci,1:trialNum)), mean(Dir_Matrix(ci,trialNum+1:trialNum*2))];
    dirIdx = find(dirM == max(dirM));
    for ss = 1:SfSzNum
        for j3 = 1:oln %%Orientaion
            rspTemporal = squeeze(G4_RspAvgOFFListTotal(ss,ci,j3,1+(dirIdx-1)*trialNum:trialNum+(dirIdx-1)*trialNum));
            rspMean = mean(rspTemporal);
            rspStd = std(rspTemporal);
            rspTrial = length(rspTemporal);
            rspMean = mean(rspTemporal);
            rspStd = std(rspTemporal);
            rspSe = rspStd./sqrt(rspTrial);
            G4_RspMeanTrialStdSeListTotal(ss,ci,j3,1) = rspMean;
            G4_RspMeanTrialStdSeListTotal(ss,ci,j3,2) = rspTrial;
            G4_RspMeanTrialStdSeListTotal(ss,ci,j3,3) = rspStd;
            G4_RspMeanTrialStdSeListTotal(ss,ci,j3,4) = rspSe;
        end
    end
end

save G4_RspMeanTrialStdSeListTotal.mat G4_RspMeanTrialStdSeListTotal;


%% Find the Maximum Response in the matrix
%%%Pick the maximum response Size in each Sf*Ori condition(3to1)

G4_SfOri_MeanSeListTotal = zeros(SfSzNum/3,cellNum,oln,2);
%1.Mean; 2.Se
G4_SzLocListTotal = zeros(SfSzNum/3,cellNum,oln);


for ci = 1:cellNum
    SfSzNum_Oln_Mean = squeeze(G4_RspMeanTrialStdSeListTotal(:,ci,:,1));
    SfSzNum_Oln_Se = squeeze(G4_RspMeanTrialStdSeListTotal(:,ci,:,4));
    for j3 = 1:oln
        for j4 = 1:SfSzNum/3
            [MaxSzRsp,MaxSz] = max(squeeze(SfSzNum_Oln_Mean(j4*3-2:j4*3,j3)));
            MaxSzSe = SfSzNum_Oln_Se(j4*3-3+MaxSz,j3);
            G4_SfOri_MeanSeListTotal(j4,ci,j3,1) = MaxSzRsp;
            G4_SfOri_MeanSeListTotal(j4,ci,j3,2) = MaxSzSe;
            G4_SzLocListTotal(j4,ci,j3) = j4*3-3+MaxSz;
        end
    end
end


save G4_SfOri_MeanSeListTotal.mat G4_SfOri_MeanSeListTotal;
save G4_SzLocListTotal.mat G4_SzLocListTotal;

%%%Pick the maximum Ori&Sf response

load G4_SfOri_MeanSeListTotal.mat;
load G4_SzLocListTotal.mat

G4_RawDataforFITTINGListTotal = zeros(cellNum,18,2);
G4_PeakSfLocListTotal = zeros(cellNum,20);
%%The 1:12 is different Ori with Max Sf, 13:18 is different Sf with Max Ori


for ci = 1:cellNum
    SfOri_Mean = squeeze(G4_SfOri_MeanSeListTotal(:,ci,:,1));
    SfOri_Se = squeeze(G4_SfOri_MeanSeListTotal(:,ci,:,2));
    SzLoc = squeeze(G4_SzLocListTotal(:,ci,:));
    [PeakSf PeakOri]=find(SfOri_Mean==max(max(SfOri_Mean)));
    G4_RawDataforFITTINGListTotal(ci,1:12,1) = SfOri_Mean(PeakSf,:);
    G4_RawDataforFITTINGListTotal(ci,13:18,1) = (SfOri_Mean(:,PeakOri))';
    G4_RawDataforFITTINGListTotal(ci,1:12,2) = SfOri_Se(PeakSf,:);
    G4_RawDataforFITTINGListTotal(ci,13:18,2) = (SfOri_Se(:,PeakOri))';
    G4_PeakSfLocListTotal(ci,1:12) = SzLoc(PeakSf,:);
    G4_PeakSfLocListTotal(ci,13) = PeakOri; %%13 is peak orientation
    G4_PeakSfLocListTotal(ci,14:19) = SzLoc(:,PeakOri);
    G4_PeakSfLocListTotal(ci,20) = PeakSf;
end


save G4_RawDataforFITTINGListTotal.mat G4_RawDataforFITTINGListTotal;
save G4_PeakSfLocListTotal.mat G4_PeakSfLocListTotal;


%% Find Anova Data

load G4_RspAvgOFFListTotal.mat;
load G4_PeakSfLocListTotal.mat;
load Y1_AnovaListTotal.mat;
load F4_AnovaListTotal.mat;
cellNum = size(G4_RspAvgOFFListTotal,2);
oln = size(G4_RspAvgOFFListTotal,3);
sfn = 6;
ol = -45:15:120;
trialNum = 6;

for k = 1:length(withdir_index)
    ci = withdir_index(k);
    dirM = [mean(Dir_Matrix(ci,1:trialNum)), mean(Dir_Matrix(ci,trialNum+1:trialNum*2))];
    dirIdx = find(dirM == max(dirM));
    anovaList = zeros(oln+sfn,trialNum);
    for j3 = 1:oln
        Loc = G4_PeakSfLocListTotal(ci,j3);
        anovaList(j3,:)= G4_RspAvgOFFListTotal(Loc,ci,j3,1+(dirIdx-1)*trialNum:trialNum+(dirIdx-1)*trialNum);
    end
    for j = oln+1:oln+sfn
        Loc2 = G4_PeakSfLocListTotal(ci,j+1);
        PeakOri = G4_PeakSfLocListTotal(ci,13);
        anovaList(j,:)= G4_RspAvgOFFListTotal(Loc2,ci,PeakOri,1+(dirIdx-1)*trialNum:trialNum+(dirIdx-1)*trialNum);
    end    
        
    anovaList = anovaList';
    anova_test = kruskalwallis(anovaList(:,1:12),1:12,'off');
    p = anova_test(1);
    Y1_AnovaListTotal(1,ci) = p; %%%OriAnova
    
    anova_test2 = kruskalwallis(anovaList(:,13:18),1:6,'off');
    p2 = anova_test2(1);
    F4_AnovaListTotal(1,ci) = p2; %%%SFAnova
end

save Y1_AnovaListTotal.mat Y1_AnovaListTotal;
save F4_AnovaListTotal.mat F4_AnovaListTotal;

% %%
% load pvalueMat.mat;
% load Y1_AnovaListTotal.mat;
% load G4_PeakSfLocListTotal.mat;
% load G4_RawDataforFITTINGListTotal.mat;
% 
% targetCellDir = intersect(find(pvalueMat<0.05), find(Y1_AnovaListTotal<0.05));
% save targetCellDir.mat targetCellDir;
% length(intersect(find(pvalueMat<0.05), find(Y1_AnovaListTotal<0.05)))/length(find(Y1_AnovaListTotal<0.05))
% targetCell = find(Y1_AnovaListTotal<0.05);
% lowIndex = [];
% targetCellOriMax = [];
% 
% for j = 1:length(targetCell)
%     targetCellOriMax = [targetCellOriMax, max(G4_RawDataforFITTINGListTotal(targetCell(j), 1:12, 1))];
%     if max(G4_RawDataforFITTINGListTotal(targetCell(j), 1:12, 1)) < 0.3
%         lowIndex = [lowIndex, j];
%     end
% end
% targetCell(lowIndex) = [];
% save targetCell.mat targetCell;
