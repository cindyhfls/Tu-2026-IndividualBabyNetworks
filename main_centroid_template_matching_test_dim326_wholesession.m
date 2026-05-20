clear
[master_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(master_dir)) % includes Util/ (smartload) and CIFTI_read_save/ (ft_read_cifti_mod)

% Load IM (network assignments)
load('IM_Tu326_19Networks.mat','IM');

% Load parcels
Parcels = smartload('Parcels_Tu_326.mat');

load('MNI_coord_meshes_32k.mat')
Anat.CtxL = MNIl;Anat.CtxR = MNIr;clear MNIl MNIr

parcelname = 'Tu_326'
load('eLABE_Y2_Y3_template_dim326_centroids.mat');
C_bestks = who('C_bestk*')
bestks = reshape(cellfun(@str2double,regexp(C_bestks,'\d+','match')),1,[]);

centroids_init_corr  = single(centroids_init_corr);
for bestk = bestks
    eval(['C_bestk',num2str(bestk),'= single(C_bestk',num2str(bestk),');'])
end
%% Assign each vertex
namestr = 'BCP_Jan2023_QCpass_asleep_atleast8min_UNC_UMN_20240124'
result_file_path = ['./template_matching_results/',namestr,'_wholesession.mat'];
T = readtable([master_dir,'/subject_tables/',namestr,'_vars.csv']);
YearGroup = T.ses_id{1}
to_overwrite = 0;

subs =  importdata([namestr,'.txt']);
if contains(namestr,'BCP')
    cohortfile = ['./cohortfiles/cohortfiles_',namestr,'_2.55sigma.txt'];
    tmasklist=['./tmasklist/tmasklist_',namestr,'_2.55sigma_QC_and_FDpt2_removeoutlierwholebrain_outliercalculatedonlowFDframes.txt'];
else
    cohortfile = ['./cohortfiles/cohortfiles_',namestr,'.txt'];
    tmasklist=['./tmasklist/tmasklist_',namestr,'.txt'];
end
% Read in subject names, functional volume locations, and surface directory
fid = fopen(cohortfile); C = textscan(fid,'%s%s%s%s'); fclose(fid);
subjects = C{1}; cifti_files = C{2};
% Read in tmasks
fid = fopen(tmasklist); C = textscan(fid,'%s%s'); fclose(fid);
tmasksubjects = C{1}; tmaskfiles = C{2};
assert(isequal(tmasksubjects,subjects),'tmasklist subjects do not match cohortfile subjects')
Nsubs = length(subs);
if exist(result_file_path,'file')
    disp('Existing results loaded');
    load(result_file_path); % so we only write the new data
    if (to_overwrite==1)
        warning('Are you sure you want to overwrite the data?')
    end
else
end
if ~exist(result_file_path,'file') || (to_overwrite ==1)
    tic
    for i = 1:Nsubs 
        fprintf('loading %i out of %i sessions...\n',i,Nsubs);
        subjectname = subs{i};
        tmask = importdata(tmaskfiles{i});
        
        if contains(T.dataset{1},'eLABE')
            [~,dtseriesname] = fileparts(cifti_files{i});
            dtseriescifti = ft_read_cifti_mod(['./datasets/eLABE/dtseries/imputedbyneighbors/',YearGroup,'/',dtseriesname,'.nii']);
        else
            dtseriescifti = ft_read_cifti_mod(cifti_files{i});
        end
        brainstructure = dtseriescifti.brainstructure(dtseriescifti.brainstructure>0);
        dtseries = single(dtseriescifti.data);
        clear dtseriescifti
        dtseries = dtseries(brainstructure<3,tmask==1);
        
        if contains(T.dataset{1},'BCP')
            ptseries_filename = dir(['./datasets/BCP/January2023/ptseries/',parcelname,'/',subjectname,'*.ptseries.nii']);
        elseif contains(T.dataset{1},'eLABE')
            ptseries_filename = dir(['./datasets/eLABE/ptseries/',YearGroup,'/',parcelname,'/',subjectname,'*.ptseries.nii']);
        end
        ptseries = ft_read_cifti_mod(fullfile(ptseries_filename.folder,ptseries_filename.name));
        ptseries = single(ptseries.data(:,tmask==1));
        
        data =  corr(dtseries',ptseries');
        
        for bestk = bestks
            eval(['C_bestk = C_bestk',num2str(bestk),';'])
            [D_top2,assn_top2] = pdist2(C_bestk,data,'correlation','Smallest',2);
            eval(['assn_all.bestk',num2str(bestk),'(:,i,:) = [assn_top2]'';']);
            eval(['D_all.bestk',num2str(bestk),'(:,i,:) = [D_top2]'';']);
        end

        [D_top2,assn_top2] = pdist2(centroids_init_corr,data,'correlation','Smallest',2); % Find the nearest centroid
        assn_all.init_corr(:,i,:) =  [assn_top2]';
        D_all.init_corr(:,i,:) = [D_top2]';
        toc
    end
    save(result_file_path,'assn_all*','D_all*');
end
return

