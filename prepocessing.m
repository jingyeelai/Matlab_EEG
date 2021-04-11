%% purpose:
%- load raw data
%- prepocess
%- store the result as variables
ccc;
%% set some parameters here
MAINPATH = '/Users/LJY/Documents/Github/';
PATHIN = [MAINPATH,'data/'];
PATHOUT = [MAINPATH,'var/'];
mkdir(PATHOUT);
cd(MAINPATH);
SUBJ = {'s1','s2','s3','s4'};                  % subject
COND = {'-with','_without','_aperiodic'};    % conditions of experiment
EVENTS={'S 10','S 20','S 30','S 40'};        % time-locking events
EP_form = -.2;                               % epoch beginning
EP_to   =  .8;                               % epoch end
BASE_from=-200;                              % baseline begin
AMP     = 100;                               % rejection threshold
hp = 1;                                      % high pass filter
hpn = 1;                                     % high pass filter after ICA
lpn = 1;                                     % low pass filter after ICA
REJ_ICA =2;                                  % pruning for ICA.
SPEED_UP= 1;                                 % 1=run ICA with PCA 20, 0= no ICA

[ALLEEG,EEG,CURRENTSET,ALLCOM] = eeglab;
%% main loop
for s = 1:length(SUBJ)
    for c = 1:length(COND)
        %  import raw data
        EEG = pop_fileio([PATHIN,SUBJ{s},'.vhdr']);
        
        % sadd channel location
        EEG = pop_chanedit(EEG, 'lookup',[MAINPATH,...
            'eeglab14_1_2b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp']);
        EEG.setname = [SUBJ{s},COND{c}];
        
        originalEEG = EEG;
        
        % re-reference to average
        EEG = pop_reref( EEG, []);
        
        % remove bad chan
        TEM = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,'ChannelCriterion',0.8,...
            'LineNoiseCriterion',4,'Highpass',hp,'BurstCriterion',20,'WindowCriterion',0.25,...
            'BurstRejection','off','Distance','Euclidian','WindowCriterionTolerances',[-Inf 7] );
        
        % interpolaration
        TEM = pop_interp(TEM, originalEEG.chanlocs, 'spherical');
        
        % epoch data into consecutive 1 sec segments
        TMP = eeg_regepochs(TMP);
        
        % reject epochs with atypical artifacts
        TMP = pop_jointprob(TMP,1,[1:length(EEG.chanlocs)] ,REJ_ICA,REJ_ICA,0,1,0,[],0);
        
        % decompose data with ICA
        if SPEED_UP
            TMP = pop_runica(TMP, 'pca',10,'interupt','on');
        else
            TMP = pop_runica(TMP, 'extended',1,'interupt','on');
        end
        
        % replace
        EEG.icawinv     = TMP.icawinv;
        EEG.icasphere   = TMP.icasphere;
        EEG.icaweights  = TMP.icaweights;
        EEG.icachansind = TMP.icachansind;
        
        % clear TMP;
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);  % ALLEEG 1
        
        % identify components representing artifacts
        pop_topoplot(EEG,0, [1:size(EEG.icawinv,2)] , EEG.setname,[6 6] ,0,'electrodes','on'); % ceil(sqrt(length(EEG.ica))) --> replace [6 6]
        pop_eegplot( EEG, 0, 1, 1);
        EEG.badcomps = input('Enter bad component indices:(if serval num. type within []):'); % use the help with IC label
        
        % back-project all other components
        EEG = pop_subcomp( EEG, EEG.badcomps, 0);
        EEG.setname = [EEG.setname,'_icacorrected' ];
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);  % ALLEEG 2
        
        % HFP & LPF
        EEG = pop_eegfiltnew(EEG, [],0.3,2750,1,[],0);
        EEG = pop_eegfiltnew(EEG, [],40,84,0,[],0);
        
        % epoch
        EEG = pop_epoch( EEG, EVENTS, [EP_form    EP_to]);
        
        % remove baseline
        EEG = pop_rmbase( EEG, [BASE_from    0]);
        EEG.setname = [SUBJ{s},COND{c} '_epoched'];
        
        % Reject epochs with residual artifacts not accounted by ICA
        EEG = pop_jointprob(EEG,1,[1:length(EEG.chanlocs)] ,3,3,0,1,0,[],0);
        EEG.setname = [EEG.setname, '_rejected'];
        EEG = pop_saveset( EEG, 'filename',[EEG.setname,'.set'],'filepath',PATHOUT);
        [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, EEG);  % ALLEEG 3
        idx = length(ALLEEG);
        
        % step 16:
        for e=1:length(EVENTS)
            EEG = pop_selectevent( ALLEEG(idx), 'latency','-2<=2','type',{EVENTS{e}},...
                'deleteevents','off','deleteepochs','on','invertepochs','off');
            EEG.setname = [EEG.setname, '_' ,CONDS{e}];
            EEG = pop_saveset( EEG, 'filename',[EEG.setname,'.set'],'filepath',PATHOUT);
            EEG = eeg_checkset(EEG);
            [ALLEEG,EEG,CURRENTSET] = eeg_store(ALLEEG, EEG);
        end
        
        % get ERP of each dataset
        ERP(s,c,:,:)=mean(EEG.data,3);
        
        
        cd(PATHOUT)
        save(['erp','s',num2str(s),'c',num2str(c),'.mat'],'ERP')
        
    end
end
eeglab redraw;