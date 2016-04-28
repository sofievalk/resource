%% ReSource Matlab Cortical Thickness DBM
% January 2016
%
% 
%
% author: bernhardt@mpg.cbs.de, valk@cbs.mpg,de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% settings  
% -----------------
% paths, toolboxes, colormaps
% -----------------
close all
clear all;

for load_data = 1

    P           = ['/scr/elbe10/boris/Resting_sofie/matlab_sofie/Project3/Boris'];
    %P           = '/Users/boris/Downloads/Boris'
    addpath([P '/toolboxes/2012-08-14_BCT/'])
    addpath([P '/toolboxes/useful/'])
    addpath([P '/toolboxes/surfstat_chicago'])
    addpath([P '/toolboxes/gretna/'])
    addpath([P '/toolboxes/NIFTI_20110921/'])
    addpath([P '/toolboxes/boris_noelsoft/boris'])
    addpath([P '/toolboxes/matlab-tools'])
    addpath([P '/toolboxes/Functions'])
    
    ico             = '5';
    kernel          = '20';
    RPATHO          = [P '/results_' kernel '/'];
    SPATH           = [P '/toolboxes/surf_fsa5'];

    mkdir(RPATHO)
    
    
    
    
    blackblue       = [zeros(1,3)*0.8;
        zeros(127,1)   (0:126)'/127   ones(127,1)];
    mycol.blackblue =  flipud(blackblue);
    
    blue            = [ones(1,3)*0.8;
        zeros(127,1)   (0:126)'/127   ones(127,1)];
    mycol.blue      =  flipud(blue);
    
    mycol.red       = [ones(1,3)*0.8; ...
        ones(64,1) linspace(0,253,64)'/254 zeros(64,1);...
        ones(64,1) ones(64,1) linspace(0,253,64)'/254];
    
    buckner         = ([102 51 51; 51 51 102; 102 102 153; 153 153 204; 153 255 51; ...
        102 204 51; 255 255 51; 255 153 51; 255 102 51; 204 51 0])./255;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %% read surface data
    % -----------------
    % mesh, parcels,
    % -----------------
    
    % load:  mesh
    % -----------------
    cd(SPATH)
    SPHERE      = SurfStatReadSurf({'lh.sphere','rh.sphere'});
    S           = SurfStatReadSurf({'lh.inflated','rh.inflated'});
    SW          = SurfStatReadSurf({'lh.white','rh.white'});
    SP          = SurfStatReadSurf({'lh.pial','rh.pial'});
    SM.tri      = SW.tri;
    SM.coord    = (SW.coord + SP.coord)./2;
    SInf        = SW;ls
    
    SInf.coord  = 0.2 *SW.coord + 0.8* S.coord;
    
    % load:  Destrieux atlas
    % -----------------
    [vertices, label, colortable] = ...
        fs_read_annotation([SPATH 'lh.aparc.a2009s.annot']);
    aparcleft = label;
    for i = 1:size(colortable.table,1)
        mycode = colortable.table(i,5);
        aparcleft(find(aparcleft == mycode)) = i;
    end
    
    [vertices, label, colortable] = ...
        fs_read_annotation([SPATH 'rh.aparc.a2009s.annot']);
    aparcright = label;
    for i = 1:size(colortable.table,1)
        mycode = colortable.table(i,5);
        aparcright(find(aparcright == mycode)) = i;
    end
    
    aparc = [aparcleft;aparcright+100];
    aparc = aparc';
    
    
    % load:  curvature
    % -----------------
    CURV = SurfStatReadData({'fsaverage_curv_lh.asc','fsaverage_curv_rh.asc'});
    
end
    
%% read csv
% -----------------
for readcsv = 1
        
        fid     = fopen([P '/time_points_scan_NOV2015_4.csv']);
        C       = textscan(fid,'%s%s%n%s%s%s%s%s%n%n%n%n%n','Delimiter',',','CollectOutput',1,'Headerlines',1);
        fclose(fid);
        dropout             = C{4}(:,4);
        bad                 = C{4}(:,5);
        tphelp              = C{3}(:,5);
        keep                = mintersect(find(~dropout), ...
            find(~bad), ...
            find(~strcmp(tphelp,'-666')), ...
            find(~strcmp(tphelp,'T4')));
        code                = C{1}(keep,1);
        age                 = C{2}(keep,1);
        sex                 = C{3}(keep,1);
        site                = C{3}(keep,2);
        group               = C{3}(keep,3);
        id                  = C{3}(keep,4);
        tp                  = C{3}(keep,5);
        year                = C{4}(keep,1);
        month               = C{4}(keep,2);
        day                 = C{4}(keep,3);
        days                = datenum(year,month,day);
        dayssince_first     = days - min(days(days>0));
        monthssince_first   = dayssince_first/30;
        
        group4              = group;
        group4(strcmp(group4,'Retest_Control'))         = {'Control'};
        group4(strcmp(group4,'Active_Control_Abends'))  = {'Group3'};
        group4(strcmp(group4,'Active_Control_Morgens')) = {'Group3'};
        
        idnum = zeros(size(id));
        tpnum = idnum;
        for i = 1:length(id)
            
            idnum(i) = str2num(id{i});
            if strcmp(tp{i},'T0')
                tpnum(i) = 0;
            elseif strcmp(tp{i},'T1')
                tpnum(i) = 1;
            elseif strcmp(tp{i},'T2')
                tpnum(i) = 2;
            elseif strcmp(tp{i},'T3')
                tpnum(i) = 3;
            end
            
        end

 
    
    %% retreat data (Appendix A) (to calculate not days between measurement but days between training
    % TC1:
    TC1_R1                = datenum(2013,6,21);
    TC1_R2                = datenum(2013,9,20);
    TC1_R3                = datenum(2014,1,3);

    % TC2: 
    TC2_R1                = datenum(2013,8,9);
    TC2_R2                = datenum(2013,11,8);
    TC2_R3                = datenum(2014,2,21);

    
    
    cd([P '/data/'])
    load('CT.mat')
    load('isthere-ct.mat')   
    load('isthere-dbm.mat')
    load('myDBMsurf.mat')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% generate change brain
    % ---------------
    zmap = load('zmyDBMsurf.mat')
    
end


for change_brain = 1
    
        minit        = zeros(size(id));
        tp_t0_there  = minit;
        tp_t1_there  = minit;
        tp_t2_there  = minit;
        tp_t3_there  = minit;
        daysfrombl   = minit;
        daysfromt1   = minit;
        daysfromt2   = minit;
        daysfromr1   = minit;
        daysfromr2   = minit;
        daysfromr3   = minit;
        CTChange_bl  = zeros(size(CT))-666;
        CTChange_t1  = zeros(size(CT))-666;
        CTChange_t2  = zeros(size(CT))-666;
        CTChange_blN = zeros(size(CT))-666;
        CTChange_T1N = zeros(size(CT))-666;
        CTChange_T2N = zeros(size(CT))-666;
        CTbl = zeros(size(CT))-666;

        DBMchange_bl = zeros(size(zmap.zdbm_map))-666;
        DBMchange_t1 = zeros(size(zmap.zdbm_map))-666;
        DBMchange_t2 = zeros(size(zmap.zdbm_map))-666;

        for i = 1:length(id)
            %% TC1
            this = find(strcmp(group{i},'Group_1'));
            if ~isempty(this)
                inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
                inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
                inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
                inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
                
                tp_t0_there(i,1)    = ~isempty(inter0) .* isthere_ct(i)==1;
                if ~isempty(inter0)
                    tp_t0_there(i,1) =  tp_t0_there(i,1) .* isthere_ct(inter0);
                end
                tp_t1_there(i,1)    = ~isempty(inter1) .* isthere_ct(i)==1;
                if ~isempty(inter1)
                    tp_t1_there(i,1) =  tp_t1_there(i,1) .* isthere_ct(inter1);
                end
                tp_t2_there(i,1)    = ~isempty(inter2) .* isthere_ct(i)==1;
                if ~isempty(inter2)
                    tp_t2_there(i,1) =  tp_t2_there(i,1) .* isthere_ct(inter2);
                end
                
                tp_t3_there(i,1)    = ~isempty(inter3) .* isthere_ct(i)==1;
                if ~isempty(inter3)
                    tp_t3_there(i,1) =  tp_t3_there(i,1) .* isthere_ct(inter3);
                end
                
                try
                    CTbl(i,:) = CT(inter0,:);
                    daysfrombl(i,1)    = days(i)-days(inter0);
                    daysfromr1(i,1)    = days(i)- TC1_R1;
                    CTChange_bl(i,:)   = CT(i,:) - CT(inter0,:);
                    CTChange_blN(i,:)  = (CT(i,:) - CT(inter0,:))./CT(inter0,:);
                    DBMchange_bl(i,:)  = zmap.zdbm_map(i,:) - zmap.zdbm_map(inter0,:);
                catch
                    %
                end
                
                try
                    daysfromt1(i,1)    = days(i)-days(inter1);
                    daysfromr2(i,1)    = days(i)-TC1_R2;
                    CTChange_t1(i,:)   = CT(i,:) - CT(inter1,:);
                    CTChange_t1N(i,:)  = (CT(i,:) - CT(inter1,:))./CT(inter1,:);
                    DBMchange_t1(i,:)  = zmap.zdbm_map(i,:) - zmap.zdbm_map(inter1,:);
                catch
                    %
                end
                
                try
                    daysfromt2(i,1)    = days(i)-days(inter2);
                    daysfromr3(i,1)    = days(i)-TC1_R3;
                    CTChange_t2(i,:)   = CT(i,:) - CT(inter2,:);
                    CTChange_t2N(i,:)  = (CT(i,:) - CT(inter2,:))./CT(inter2,:);
                    DBMchange_t2(i,:)  = zmap.zdbm_map(i,:) - zmap.zdbm_map(inter2,:);
                    
                catch
                    %
                end
            else
                inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
                inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
                inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
                inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
                
                tp_t0_there(i,1)    = ~isempty(inter0) .* isthere_ct(i)==1;
                if ~isempty(inter0)
                    tp_t0_there(i,1) =  tp_t0_there(i,1) .* isthere_ct(inter0);
                end
                tp_t1_there(i,1)    = ~isempty(inter1) .* isthere_ct(i)==1;
                if ~isempty(inter1)
                    tp_t1_there(i,1) =  tp_t1_there(i,1) .* isthere_ct(inter1);
                end
                tp_t2_there(i,1)    = ~isempty(inter2) .* isthere_ct(i)==1;
                if ~isempty(inter2)
                    tp_t2_there(i,1) =  tp_t2_there(i,1) .* isthere_ct(inter2);
                end
                
                tp_t3_there(i,1)    = ~isempty(inter3) .* isthere_ct(i)==1;
                if ~isempty(inter3)
                    tp_t3_there(i,1) =  tp_t3_there(i,1) .* isthere_ct(inter3);
                end
                
                try
                    CTbl(i,:) = CT(inter0,:);
                    daysfrombl(i,1)    = days(i)-days(inter0);
                    daysfromr1(i,1)    = days(i)- TC2_R1;
                    CTChange_bl(i,:)   = CT(i,:) - CT(inter0,:);
                    CTChange_blN(i,:)  = (CT(i,:) - CT(inter0,:))./CT(inter0,:);
                    DBMchange_bl(i,:)  = zmap.zdbm_map(i,:) - zmap.zdbm_map(inter0,:);
                catch
                    %
                end
                
                try
                    daysfromt1(i,1)    = days(i)-days(inter1);
                    daysfromr2(i,1)    = days(i)-TC2_R2;
                    CTChange_t1(i,:)   = CT(i,:) - CT(inter1,:);
                    CTChange_t1N(i,:)  = (CT(i,:) - CT(inter1,:))./CT(inter1,:);
                    DBMchange_t1(i,:)  = zmap.zdbm_map(i,:) - zmap.zdbm_map(inter1,:);
                catch
                    %
                end
                
                try
                    daysfromt2(i,1)    = days(i)-days(inter2);
                    daysfromr3(i,1)    = days(i)-TC2_R3;
                    CTChange_t2(i,:)   = CT(i,:) - CT(inter2,:);
                    CTChange_t2N(i,:)  = (CT(i,:) - CT(inter2,:))./CT(inter2,:);
                    DBMchange_t2(i,:)  = zmap.zdbm_map(i,:) - zmap.zdbm_map(inter2,:);
                    
                catch
                    %
                end
                
            end
        end
                
                
                %     dat = [ idnum/10000 tpnum isthere_ct tp_t0_there tp_t1_there tp_t2_there ...
                %         daysfrombl./max(abs(daysfrombl)) daysfromt1./max(abs(daysfromt1)) ...
                %         daysfromt2./max(abs(daysfromt2))];
                
                % f=figure, imagesc(dat)
                mask = mean(CT,1)>0.9;
        
end
  
for load_behav = 1 
    
    cd([P '/data/'])
    load('P_body.mat')
    load('P_breath.mat')
    load('A_hearth.mat')
    load('A_adyad.mat')
    load('Pe_thought.mat')
    load('Pe_pdyad.mat')
    
    
    load('tom.mat')
    load('compassion.mat')
    load('attention_c.mat')
    load('attention.mat')
    
    

    % not all subjects are there for the resting!!
    load('rp_arest.mat')
    load('FD_Power_RS.mat')
    
end

for initi = 1
    
    minit               = zeros(size(id))-666;
    
    com_bl              = minit;
    com_bl2             = minit;
    com_bl3             = minit;
    comChange_bl        = minit;
    comChange_t1        = minit;
    comChange_t2        = minit;
    
    comChange_blN        = minit;
    comChange_t1N        = minit;
    comChange_t2N        = minit;
    
    att_bl              = minit;
    att_bl2              = minit;
    att_bl3              = minit;

    attChange_bl        = minit;
    attChange_t1        = minit;
    attChange_t2        = minit;
    
    attChange_blN        = minit;
    attChange_t1N        = minit;
    attChange_t2N        = minit;
    
    att_bl_c              = minit;
    attChange_bl_c        = minit;
    attChange_t1_c        = minit;
    attChange_t2_c        = minit;
    
    tom_bl              = minit;
    tom_bl2              = minit;
    tom_bl3              = minit;
    tomChange_bl        = minit;
    tomChange_t1        = minit;
    tomChange_t2        = minit;
    
    tomChange_blN        = minit;
    tomChange_t1N        = minit;
    tomChange_t2N        = minit;
    
    hearChange_bl        = minit;
    hearChange_t1        = minit;
    hearChange_t2        = minit;

    rp_arestChange_bl    = zeros(size(FD_Power_RS))-666;
    rp_arestChange_t1    = zeros(size(FD_Power_RS))-666;
    rp_arestChange_t2    = zeros(size(FD_Power_RS))-666;
    
end
    
    
%% attention
for i = 1: length(id)
    i
    inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
    inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
    inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
    inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
    
    attention2 = (attention);
    if ~isempty(inter1) & attention(inter1)> -500
        if ~isempty(inter0) & attention(inter0)> -500
            att_bl(i)        = (attention(inter0,1));
            
            attChange_bl(i)  = (attention2(i)-attention2(inter0));
            attChange_blN(i) = (attention(i)-attention(inter0))./attention(inter0);
            attChange_bf(i)  =  attention(inter0);
            
            attChange_bl_c(i)  = (attention_c(i)-attention_c(inter0));
            attChange_blN_c(i) = (attention_c(i)-attention_c(inter0))./attention_c(inter0);
          
        end
    end
    
    if ~isempty(inter2) & attention(inter2)> -500
        if ~isempty(inter1)  & attention(inter1)> -500
            att_bl2(i)        = (attention(inter1,1));
            attChange_t1(i)  = (attention2(i)-attention2(inter1));
            attChange_t1N(i) = (attention(i)-attention(inter1))./attention(inter1);
            
            attChange_t1_c(i)  = (attention_c(i)-attention_c(inter1));
            attChange_t1N_c(i) = (attention_c(i)-attention_c(inter1))./attention_c(inter1);
            
          
            
        end
    end
    
    if ~isempty(inter3)  & attention(inter3)> -500
        if ~isempty(inter2)  & attention(inter2)> -500
            att_bl3(i)        = (attention(inter2,1));
            attChange_t2(i)  = (attention2(i)-attention2(inter2));
            attChange_t2N(i) = (attention(i)-attention(inter2))./attention(inter2);
            
            attChange_t2_c(i)  = (attention_c(i)-attention_c(inter2));
            attChange_t2N_c(i) = (attention_c(i)-attention_c(inter2))./attention_c(inter2);
            
          
        end
    end
end

%% tom
for i = 1: length(id)
    i
    inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
    inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
    inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
    inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
    
    
    if ~isempty(inter1) & tom(inter1)> -500
        if ~isempty(inter0) & tom(inter0)> -500
            tom_bl(i)    = tom(inter0);
            tomChange_bl(i)  = (tom(i)-tom(inter0));
            tomChange_blN(i) = (tom(i)-tom(inter0))./tom(inter0);
            
            hearChange_bl(i)  = (hearing(i)-hearing(inter0));
        end
    end
    
    if ~isempty(inter2) & tom(inter2)> -500
        if ~isempty(inter1) & tom(inter1)> -500
            tom_bl2(i)    = tom(inter1);
            tomChange_t1(i)  = (tom(i)-tom(inter1));
            tomChange_t1N(i) = (tom(i)-tom(inter1))./tom(inter1);
            tomChange_t1b(i) =  tom(inter1);
            
            hearChange_t1(i)  = (hearing(i)-hearing(inter1));
            
        end
    end
    
    if ~isempty(inter3) & tom(inter3)> -500
        if ~isempty(inter2) & tom(inter2)> -500
            tom_bl3(i)    = tom(inter2);
            tomChange_t2(i)  = (tom(i)-tom(inter2));
            tomChange_t2N(i) = (tom(i)-tom(inter2))./tom(inter2);
            
            hearChange_t2(i)  = (hearing(i)-hearing(inter2));
        end
    end
    
end

%% compassion
for i = 1: length(id)
    i
    inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
    inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
    inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
    inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
    
    compassion2 = (compassion);
    if ~isempty(inter1) & compassion(inter1)> -500
        if ~isempty(inter0) & compassion(inter0)> -500
            com_bl(i)  = compassion(inter0);
            comChange_bl(i)  = (compassion2(i)-compassion2(inter0));
            comChange_blN(i) = (compassion(i)-compassion(inter0))./compassion(inter0);
            comChange_bf(i)  = compassion(inter0);
            
        end
    end
    
    if ~isempty(inter2) & compassion(inter2)> -500
        if ~isempty(inter1) & compassion(inter1)> -500
            com_bl2(i)  = compassion(inter1);
            comChange_t1(i)  = (compassion2(i)-compassion2(inter1));
            comChange_t1N(i) = (compassion(i)-compassion(inter1))./compassion(inter1);
            
        end
    end
    
    if ~isempty(inter3) & compassion(inter3)> -500
        if ~isempty(inter2) & compassion(inter2)> -500
            com_bl3(i)  = compassion(inter2);
            comChange_t2(i)  = (compassion2(i)-compassion2(inter2));
            comChange_t2N(i) = (compassion(i)-compassion(inter2))./compassion(inter2);
            
        end
    end
    
end

%% resting-move
for i = 1: length(id)
    i
    inter0 		  = intersect(find(tpnum==0),find(idnum==idnum(i)));
    inter1        = intersect(find(tpnum==1),find(idnum==idnum(i)));
    inter2        = intersect(find(tpnum==2),find(idnum==idnum(i)));
    inter3        = intersect(find(tpnum==3),find(idnum==idnum(i)));
    

    if ~isempty(inter1) & FD_Power_RS(inter1,1)~=0
        if ~isempty(inter0) & FD_Power_RS(inter0,1)~=0
            rp_arestChange_bl(i,:)  = (FD_Power_RS(i,:)-FD_Power_RS(inter0,:));       
        end
    end
    
    if ~isempty(inter2) & FD_Power_RS(inter2,1)~=0
        if ~isempty(inter1) & FD_Power_RS(inter1,1)~=0
            rp_arestChange_t1(i,:)  = (FD_Power_RS(i,:)-FD_Power_RS(inter1,:));                 
        end
    end
    
    if ~isempty(inter3) & FD_Power_RS(inter3,1)~=0
        if ~isempty(inter2) & FD_Power_RS(inter2,1)~=0
            rp_arestChange_t2(i,:)  = (FD_Power_RS(i,:)-FD_Power_RS(inter2,:));                 
        end
    end
    
end

%% make behavior and groups re-ordered
% ---------------
% label the training blocks
% ---------------
for reorder = 1
    groupN = group4;
    groupN(strcmp(group4,'Group_1') & strcmp(tp,'T1'))          = {'Presence'};
    groupN(strcmp(group4,'Group_1') & strcmp(tp,'T2'))          = {'Affect'};
    groupN(strcmp(group4,'Group_1') & strcmp(tp,'T3'))          = {'Perspective'};
    
    groupN(strcmp(group4,'Group_2') & strcmp(tp,'T1'))          = {'Presence'};
    groupN(strcmp(group4,'Group_2') & strcmp(tp,'T2'))          = {'Perspective'};
    groupN(strcmp(group4,'Group_2') & strcmp(tp,'T3'))          = {'Affect'};
    
    groupN(strcmp(group4,'Group3')  & strcmp(tp,'T1'))          = {'Affect'};
    
    groupN(strcmp(group4,'Control') & strcmp(tp,'T1'))          = {'Control1'};
    groupN(strcmp(group4,'Control') & strcmp(tp,'T2'))          = {'Control1'};
    groupN(strcmp(group4,'Control') & strcmp(tp,'T3'))          = {'Control1'};
    
    
    %% remodel brain change wrt last time point
    % ---------------
    % reoder differences in thickness
    % ---------------
    CTChange_last = zeros(length(id),size(CTChange_bl,2));
    CTChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:) = CTChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:);
    CTChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:) = CTChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:);
    CTChange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:) = CTChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:);
    
    CTChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:) = CTChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:);
    CTChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:) = CTChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:);
    CTChange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:) = CTChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:);
    
    CTChange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:)  = CTChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:);
    CTChange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:)  = CTChange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:);
    CTChange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:)  = CTChange_t2(strcmp(group4,'Group3') & strcmp(tp,'T3'),:);
    
    CTChange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:) = CTChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:);
    CTChange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:) = CTChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:);
    CTChange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:) = CTChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:);
    
    % for T0-T2 differences in thickness...
    CTChangeBL = zeros(length(id),size(CTChange_bl,2));
    CTChangeBL(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:) = CTChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:);
    CTChangeBL(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:) = CTChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:);
    CTChangeBL(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:) = CTChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:);
    
    CTChangeBL(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:) = CTChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:);
    CTChangeBL(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:) = CTChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:);
    CTChangeBL(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:) = CTChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:);
    
    CTChangeBL(strcmp(group4,'Group3') & strcmp(tp,'T1'),:)  = CTChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:);
    CTChangeBL(strcmp(group4,'Group3') & strcmp(tp,'T2'),:)  = CTChange_bl(strcmp(group4,'Group3') & strcmp(tp,'T2'),:);
    CTChangeBL(strcmp(group4,'Group3') & strcmp(tp,'T3'),:)  = CTChange_bl(strcmp(group4,'Group3') & strcmp(tp,'T3'),:);
    
    CTChangeBL(strcmp(group4,'Control') & strcmp(tp,'T1'),:) = CTChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:);
    CTChangeBL(strcmp(group4,'Control') & strcmp(tp,'T2'),:) = CTChange_bl(strcmp(group4,'Control') & strcmp(tp,'T2'),:);
    CTChangeBL(strcmp(group4,'Control') & strcmp(tp,'T3'),:) = CTChange_bl(strcmp(group4,'Control') & strcmp(tp,'T3'),:);
    
    
    % % normalized differences in thickness
    % % ---------------
    CTChange_lastN = zeros(length(id),size(CTChange_blN,2));
    CTChange_lastN(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:) = CTChange_blN(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:);
    CTChange_lastN(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:) = CTChange_t1N(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:);
    CTChange_lastN(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:) = CTChange_t2N(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:);
    
    CTChange_lastN(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:) = CTChange_blN(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:);
    CTChange_lastN(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:) = CTChange_t1N(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:);
    CTChange_lastN(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:) = CTChange_t2N(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:);
    
    CTChange_lastN(strcmp(group4,'Group3') & strcmp(tp,'T1'),:)  = CTChange_blN(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:);
    CTChange_lastN(strcmp(group4,'Group3') & strcmp(tp,'T2'),:)  = CTChange_t1N(strcmp(group4,'Group3') & strcmp(tp,'T2'),:);
    CTChange_lastN(strcmp(group4,'Group3') & strcmp(tp,'T3'),:)  = CTChange_t2N(strcmp(group4,'Group3') & strcmp(tp,'T3'),:);
    
    CTChange_lastN(strcmp(group4,'Control') & strcmp(tp,'T1'),:) = CTChange_blN(strcmp(group4,'Control') & strcmp(tp,'T1'),:);
    CTChange_lastN(strcmp(group4,'Control') & strcmp(tp,'T2'),:) = CTChange_t1N(strcmp(group4,'Control') & strcmp(tp,'T2'),:);
    CTChange_lastN(strcmp(group4,'Control') & strcmp(tp,'T3'),:) = CTChange_t2N(strcmp(group4,'Control') & strcmp(tp,'T3'),:);
    
    
    %% remodel dbm change wrt last time point
    DBMchange_last = zeros(length(id),size(DBMchange_bl,2));
    
    DBMchange_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:) = DBMchange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:);
    DBMchange_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:) = DBMchange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:);
    DBMchange_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:) = DBMchange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:);
    
    DBMchange_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:) = DBMchange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:);
    DBMchange_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:) = DBMchange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:);
    DBMchange_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:) = DBMchange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:);
    
    
    DBMchange_last(strcmp(group4,'Group3') & strcmp(tp,'T1'),:) = DBMchange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:);
    DBMchange_last(strcmp(group4,'Group3') & strcmp(tp,'T2'),:)  = DBMchange_t1(strcmp(group4,'Group3') & strcmp(tp,'T2'),:);
    DBMchange_last(strcmp(group4,'Group3') & strcmp(tp,'T3'),:)  = DBMchange_t1(strcmp(group4,'Group3') & strcmp(tp,'T3'),:);
    
    DBMchange_last(strcmp(group4,'Control') & strcmp(tp,'T1'),:) = DBMchange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:);
    DBMchange_last(strcmp(group4,'Control') & strcmp(tp,'T2'),:) = DBMchange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:);
    DBMchange_last(strcmp(group4,'Control') & strcmp(tp,'T3'),:) = DBMchange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:);
    
    
    %% remodel days wrt last time point
    Days_last = zeros(length(id),1);
    Days_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'))        = daysfrombl(strcmp(group4,'Group_1') & strcmp(tp,'T1'));
    Days_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'))        = daysfromt1(strcmp(group4,'Group_1') & strcmp(tp,'T2'));
    Days_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'))        = daysfromt2(strcmp(group4,'Group_1') & strcmp(tp,'T3'));
    
    Days_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'))        = daysfrombl(strcmp(group4,'Group_2') & strcmp(tp,'T1'));
    Days_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'))        = daysfromt1(strcmp(group4,'Group_2') & strcmp(tp,'T2'));
    Days_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'))        = daysfromt2(strcmp(group4,'Group_2') & strcmp(tp,'T3'));
    
    Days_last(strcmp(group4,'Control') & strcmp(tp,'T1'))        = daysfrombl(strcmp(group4,'Control') & strcmp(tp,'T1'));
    Days_last(strcmp(group4,'Control') & strcmp(tp,'T2'))        = daysfromt1(strcmp(group4,'Control') & strcmp(tp,'T2'));
    Days_last(strcmp(group4,'Control') & strcmp(tp,'T3'))        = daysfromt2(strcmp(group4,'Control') & strcmp(tp,'T3'));
    
    Days_last(strcmp(group4,'Group3')  & strcmp(tp,'T1'))        = daysfrombl(strcmp(group4,'Group3')  & strcmp(tp,'T1'));
    Days_last(strcmp(group4,'Group3') & strcmp(tp,'T2'))         = daysfromt1(strcmp(group4,'Group3') & strcmp(tp,'T2'));
    Days_last(strcmp(group4,'Group3') & strcmp(tp,'T3'))         = daysfromt2(strcmp(group4,'Group3') & strcmp(tp,'T3'));
    
    %% remodel days wrt last time point- retreats
    Days_last = zeros(length(id),1);
    Days_last(strcmp(group4,'Group_1') & strcmp(tp,'T1'))        = daysfromr1(strcmp(group4,'Group_1') & strcmp(tp,'T1'));
    Days_last(strcmp(group4,'Group_1') & strcmp(tp,'T2'))        = daysfromr2(strcmp(group4,'Group_1') & strcmp(tp,'T2'));
    Days_last(strcmp(group4,'Group_1') & strcmp(tp,'T3'))        = daysfromr3(strcmp(group4,'Group_1') & strcmp(tp,'T3'));
    
    Days_last(strcmp(group4,'Group_2') & strcmp(tp,'T1'))        = daysfromr1(strcmp(group4,'Group_2') & strcmp(tp,'T1'));
    Days_last(strcmp(group4,'Group_2') & strcmp(tp,'T2'))        = daysfromr2(strcmp(group4,'Group_2') & strcmp(tp,'T2'));
    Days_last(strcmp(group4,'Group_2') & strcmp(tp,'T3'))        = daysfromr3(strcmp(group4,'Group_2') & strcmp(tp,'T3'));
    
    Days_last(strcmp(group4,'Control') & strcmp(tp,'T1'))        = daysfrombl(strcmp(group4,'Control') & strcmp(tp,'T1'));
    Days_last(strcmp(group4,'Control') & strcmp(tp,'T2'))        = daysfromt1(strcmp(group4,'Control') & strcmp(tp,'T2'));
    Days_last(strcmp(group4,'Control') & strcmp(tp,'T3'))        = daysfromt2(strcmp(group4,'Control') & strcmp(tp,'T3'));
    
    Days_last(strcmp(group4,'Group3')  & strcmp(tp,'T1'))        = daysfrombl(strcmp(group4,'Group3')  & strcmp(tp,'T1'));
    Days_last(strcmp(group4,'Group3') & strcmp(tp,'T2'))         = daysfromt1(strcmp(group4,'Group3') & strcmp(tp,'T2'));
    Days_last(strcmp(group4,'Group3') & strcmp(tp,'T3'))         = daysfromt2(strcmp(group4,'Group3') & strcmp(tp,'T3'));
   
    %% remake behavioral change wrt last time point
    %--- attention change % note the slight difference in syntax
    chA = zeros(length(id),1);
    chA(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:)           = attChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),1);
    chA(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:)           = attChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),1);
    chA(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:)           = attChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),1);
    
    chA(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:)           = attChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),1);
    chA(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:)           = attChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),1);
    chA(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:)           = attChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),1);
    
    chA(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:)           = attChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),1);
    
    chA(strcmp(group4,'Control') & strcmp(tp,'T1'),:)           = attChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),1);
    chA(strcmp(group4,'Control') & strcmp(tp,'T2'),:)           = attChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),1);
    chA(strcmp(group4,'Control') & strcmp(tp,'T3'),:)           = attChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),1);
    
    chA_c = zeros(length(id),1);
    chA_c(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:)           = attChange_bl_c(strcmp(group4,'Group_1') & strcmp(tp,'T1'),1);
    chA_c(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:)           = attChange_t1_c(strcmp(group4,'Group_1') & strcmp(tp,'T2'),1);
    chA_c(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:)           = attChange_t2_c(strcmp(group4,'Group_1') & strcmp(tp,'T3'),1);
    
    chA_c(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:)           = attChange_bl_c(strcmp(group4,'Group_2') & strcmp(tp,'T1'),1);
    chA_c(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:)           = attChange_t1_c(strcmp(group4,'Group_2') & strcmp(tp,'T2'),1);
    chA_c(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:)           = attChange_t2_c(strcmp(group4,'Group_2') & strcmp(tp,'T3'),1);
    
    chA_c(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:)           = attChange_bl_c(strcmp(group4,'Group3')  & strcmp(tp,'T1'),1);
    
    chA_c(strcmp(group4,'Control') & strcmp(tp,'T1'),:)           = attChange_bl_c(strcmp(group4,'Control') & strcmp(tp,'T1'),1);
    chA_c(strcmp(group4,'Control') & strcmp(tp,'T2'),:)           = attChange_t1_c(strcmp(group4,'Control') & strcmp(tp,'T2'),1);
    chA_c(strcmp(group4,'Control') & strcmp(tp,'T3'),:)           = attChange_t2_c(strcmp(group4,'Control') & strcmp(tp,'T3'),1);
    
    %%baseline
    chAbl = zeros(length(id),1);
    chAbl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:)           = att_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),1);
    chAbl(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:)           = att_bl2(strcmp(group4,'Group_1') & strcmp(tp,'T2'),1);
    chAbl(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:)           = att_bl3(strcmp(group4,'Group_1') & strcmp(tp,'T3'),1);
    
    chAbl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:)           = att_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),1);
    chAbl(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:)           = att_bl2(strcmp(group4,'Group_2') & strcmp(tp,'T2'),1);
    chAbl(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:)           = att_bl3(strcmp(group4,'Group_2') & strcmp(tp,'T3'),1);
    
    chAbl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:)           = att_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),1);
    
    chAbl(strcmp(group4,'Control') & strcmp(tp,'T1'),:)           = att_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),1);
    chAbl(strcmp(group4,'Control') & strcmp(tp,'T2'),:)           = att_bl2(strcmp(group4,'Control') & strcmp(tp,'T2'),1);
    chAbl(strcmp(group4,'Control') & strcmp(tp,'T3'),:)           = att_bl3(strcmp(group4,'Control') & strcmp(tp,'T3'),1);
   
    chc = zeros(length(id),1);
    chc(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:)           = comChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),1);
    chc(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:)           = comChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),1);
    chc(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:)           = comChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),1);
    
    chc(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:)           = comChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),1);
    chc(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:)           = comChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),1);
    chc(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:)           = comChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),1);
    
    chc(strcmp(group4,'Group3') &  strcmp(tp,'T1'),:)           = comChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),1);
    
    chc(strcmp(group4,'Control') & strcmp(tp,'T1'),:)           = comChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),1);
    chc(strcmp(group4,'Control') & strcmp(tp,'T2'),:)           = comChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),1);
    chc(strcmp(group4,'Control') & strcmp(tp,'T3'),:)           = comChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),1);
    
    %      %%baseline
    chcbl = zeros(length(id),1);
    chcbl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:)           = com_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),1);
    chcbl(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:)           = com_bl2(strcmp(group4,'Group_1') & strcmp(tp,'T2'),1);
    chcbl(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:)           = com_bl3(strcmp(group4,'Group_1') & strcmp(tp,'T3'),1);
    
    chcbl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:)           = com_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),1);
    chcbl(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:)           = com_bl2(strcmp(group4,'Group_2') & strcmp(tp,'T2'),1);
    chcbl(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:)           = com_bl3(strcmp(group4,'Group_2') & strcmp(tp,'T3'),1);
    
    chcbl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:)           = com_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),1);
    
    chcbl(strcmp(group4,'Control') & strcmp(tp,'T1'),:)           = com_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),1);
    chcbl(strcmp(group4,'Control') & strcmp(tp,'T2'),:)           = com_bl2(strcmp(group4,'Control') & strcmp(tp,'T2'),1);
    chcbl(strcmp(group4,'Control') & strcmp(tp,'T3'),:)           = com_bl3(strcmp(group4,'Control') & strcmp(tp,'T3'),1);
   
    %--- tom change
    cht = zeros(length(id),1);
    cht(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:)           = tomChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),1);
    cht(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:)           = tomChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),1);
    cht(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:)           = tomChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),1);
    
    cht(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:)           = tomChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),1);
    cht(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:)           = tomChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),1);
    cht(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:)           = tomChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),1);
    
    cht(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:)           = tomChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),1);
    
    cht(strcmp(group4,'Control') & strcmp(tp,'T1'),:)           = tomChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),1);
    cht(strcmp(group4,'Control') & strcmp(tp,'T2'),:)           = tomChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),1);
    cht(strcmp(group4,'Control') & strcmp(tp,'T3'),:)           = tomChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),1);
    
%      %%baseline
    chtbl = zeros(length(id),1);
    chtbl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:)           = tom_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),1);
    chtbl(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:)           = tom_bl2(strcmp(group4,'Group_1') & strcmp(tp,'T2'),1);
    chtbl(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:)           = tom_bl3(strcmp(group4,'Group_1') & strcmp(tp,'T3'),1);
    
    chtbl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:)           = tom_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),1);
    chtbl(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:)           = tom_bl2(strcmp(group4,'Group_2') & strcmp(tp,'T2'),1);
    chtbl(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:)           = tom_bl3(strcmp(group4,'Group_2') & strcmp(tp,'T3'),1);
    
    chtbl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:)           = tom_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),1);
    
    chtbl(strcmp(group4,'Control') & strcmp(tp,'T1'),:)           = tom_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),1);
    chtbl(strcmp(group4,'Control') & strcmp(tp,'T2'),:)           = tom_bl2(strcmp(group4,'Control') & strcmp(tp,'T2'),1);
    chtbl(strcmp(group4,'Control') & strcmp(tp,'T3'),:)           = tom_bl3(strcmp(group4,'Control') & strcmp(tp,'T3'),1);
   
  %     %--- movement change
    mvt = zeros(length(id),1);
    mvt(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:)           = rp_arestChange_bl(strcmp(group4,'Group_1') & strcmp(tp,'T1'),:);
    mvt(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:)           = rp_arestChange_t1(strcmp(group4,'Group_1') & strcmp(tp,'T2'),:);
    mvt(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:)           = rp_arestChange_t2(strcmp(group4,'Group_1') & strcmp(tp,'T3'),:);
    
    mvt(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:)           = rp_arestChange_bl(strcmp(group4,'Group_2') & strcmp(tp,'T1'),:);
    mvt(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:)           = rp_arestChange_t1(strcmp(group4,'Group_2') & strcmp(tp,'T2'),:);
    mvt(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:)           = rp_arestChange_t2(strcmp(group4,'Group_2') & strcmp(tp,'T3'),:);
    
    mvt(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:)           = rp_arestChange_bl(strcmp(group4,'Group3')  & strcmp(tp,'T1'),:);
    
    mvt(strcmp(group4,'Control') & strcmp(tp,'T1'),:)           = rp_arestChange_bl(strcmp(group4,'Control') & strcmp(tp,'T1'),:);
    mvt(strcmp(group4,'Control') & strcmp(tp,'T2'),:)           = rp_arestChange_t1(strcmp(group4,'Control') & strcmp(tp,'T2'),:);
    mvt(strcmp(group4,'Control') & strcmp(tp,'T3'),:)           = rp_arestChange_t2(strcmp(group4,'Control') & strcmp(tp,'T3'),:);
    
    
    %% data
    X                   = CTChange_last;
    X(isnan(X))         = 0;
    X(isinf(X))         = 0;
    
    XBL                 = CTChangeBL;
    XBL(isnan(XBL))     = 0;
    XBL(isinf(XBL))     = 0;
    
    XB                  = DBMchange_last;
    XB(isnan(XB))       = 0;
    XB(isinf(XB))       = 0;

    Att                 = chA;
    Att(abs(Att)>100)   = -666;
    
    Att_c                 = chA_c;
    Att_c(abs(Att_c)>100)   = -666;
    
    Comp                = chc;
    Comp(abs(chc)>100)  = -666;
    
    Tom                 = cht;
    Tom(abs(Tom)>100)   = -666;
    
    Move                = mvt;
    Move(abs(Move)>100) = -666;
end


%% make ct baseline;


%% //Preparestuff ---------------------------------------------------------

%% analysis settings
% -----------------
% these are the analysis settings for the Resource analysis 

p_threshs       = [0.001 0.005 0.01 0.025 0.03];
lobescheck      = [1 2 3 4 5 6 9 10];   
pval            = 4; 

p_thresh        = p_threshs(pval);    % cluster-level threshold 
RPATH           = [RPATHO '/APRIL_' num2str(p_thresh) '/'];
mkdir(RPATH)
pwe = 0.025; 


%% SUBJECTS INCLUSIONS
for subjects_incl = 1
    keep = find(isthere_ct ==1), find(~strcmp(group4,'Group3')))
    size(keep) %N=848
    size((find(tpnum(keep)==0))) % N = 232
    size(intersect(find(tpnum(keep)==0), find(strcmp(group4(keep),'Group_1')))) %77/80 
    size(intersect(find(tpnum(keep)==0), find(strcmp(group4(keep),'Group_2')))) %74/81 
    size(intersect(find(tpnum(keep)==0), find(strcmp(group4(keep),'Control')))) %81/90 
    idk = id(keep)
    missing_T0_G1 = idk(intersect(find(tpnum(keep)==0), find(strcmp(group(keep),'Group_2'))))
    missing_T0_G2 = idk(intersect(find(tpnum(keep)==0), find(strcmp(group(keep),'Group_2'))))
    missing_T0_RCC = idk(intersect(find(tpnum(keep)==0), find(strcmp(group4(keep),'Control'))))
 
    mean(age) %40.52
    std(age)  %9.3

    size((find(tpnum(keep)==1))) % N = 213
    size(intersect(find(tpnum(keep)==1), find(strcmp(group4(keep),'Group_1')))) %69/80
    size(intersect(find(tpnum(keep)==1), find(strcmp(group4(keep),'Group_2')))) %69/81
    size(intersect(find(tpnum(keep)==1), find(strcmp(group4(keep),'Control')))) %75/90

    missing_T1_G1 = idk(intersect(find(tpnum(keep)==1), find(strcmp(group(keep),'Group_1')))) 
    missing_T1_G2 = idk(intersect(find(tpnum(keep)==1), find(strcmp(group(keep),'Group_2'))));
    missing_T1_RCC = idk(intersect(find(tpnum(keep)==1), find(strcmp(group4(keep),'Control'))));

    size((find(tpnum(keep)==2))) % N = 204
    size(intersect(find(tpnum(keep)==2), find(strcmp(group4(keep),'Group_1')))) %65
    size(intersect(find(tpnum(keep)==2), find(strcmp(group4(keep),'Group_2')))) %69
    size(intersect(find(tpnum(keep)==2), find(strcmp(group4(keep),'Control')))) %70

    size((find(tpnum(keep)==3))) % N = 199
    size(intersect(find(tpnum(keep)==3), find(strcmp(group4(keep),'Group_1')))) %59
    size(intersect(find(tpnum(keep)==3), find(strcmp(group4(keep),'Group_2')))) %68
    size(intersect(find(tpnum(keep)==3), find(strcmp(group4(keep),'Control')))) %72

    
    keep    = mintersect(find(isthere_ct ==1),...
        find(tp_t0_there==1), ...    
        find(tp_t1_there==1), ...
        find(tp_t2_there==1), ...
        find(tp_t3_there==1), ...
        union(find(tpnum==0), union(find(tpnum==1),union(find(tpnum==3), find(tpnum==2 )))));
    
    size(keep) %N = 716
    size((find(tpnum(keep)==0))) % N = 179
    size(intersect(find(tpnum(keep)==0), find(strcmp(group4(keep),'Group_1')))) %54
    size(intersect(find(tpnum(keep)==0), find(strcmp(group4(keep),'Group_2')))) %59
    size(intersect(find(tpnum(keep)==0), find(strcmp(group4(keep),'Control')))) %66
    
    mean(age(keep)) %40.3
    std(age(keep))  %9.6
    size(intersect(find(tpnum(keep)==0), find(strcmp(sex(keep),'female')))) %96 (81)
end

%% for QC

for QC = 1
    istherekeep = find(isthere_ct);
    for i= 1:length(istherekeep)
        f = figure,
        BoSurfStatViewData(CT(istherekeep(i),:), SInf,'')
        exportfigbo(f,[RPATH 'sub' num2str(idnum(keep(i))) '_T' num2str(tpnum(keep(i))) 'CT.png'],'png',10);
        close(f)
    end
    
     keep    = mintersect(find(isthere_dbm), ...
        find(tp_t3_there==1), ...
        find(tp_t1_there==1), ...
        find(tp_t2_there==1 ), ...
        find(tp_t0_there==1 ), ...
        union(find(tpnum==0), union(find(tpnum==1),union(find(tpnum==3), find(tpnum==2 )))));
  
    istherekeep = keep;
    for i= 1 :length(istherekeep)
        f = figure,
        BoSurfStatViewData(zmap.zdbm_map(istherekeep(i),:), SInf,'')
        BoSurfStatColLim([-3 3])
        exportfigbo(f,[RPATH 'sub' num2str(idnum(istherekeep(i))) '_T' num2str(tpnum(istherekeep(i))) 'DBM.png'],'png',10);
        close(f)
    end
end

%% load AAL atlas


for aal = 1
    cd([P '/data/')
    load('aal.mat')
end


%% CT CHANGE FIGURE 2 and 3 | SUPPLEMENTs
for FIG12 = 1
    keep    = mintersect(find(isthere_ct), ...
        find(tp_t0_there==1), ...
        find(tp_t3_there==1), ...
        find(tp_t1_there==1), ...
        find(tp_t2_there==1 ), ...
        union(find(tpnum==0), union(find(tpnum==1),union(find(tpnum==3), find(tpnum==2 )))));
    keepMAIN = keep;
    size(keep) % N = 716
    mask    = mean(CT,1)>0.4;
    Ck      = CT(keep,:);
    ik      = id(keep,:);
    ink     = idnum(keep);
    tk      = tpnum(keep,:);
    dk      = Days_last(keep);
    dkfac   = tp(keep);
    gk      = group(keep);
    g4k     = group4(keep);
    ak      = age(keep);
    sk      = sex(keep);
    
    A       = term(ak);
    S       = term(sk);
    G       = term(gk);
    G4      = term(g4k);
    D       = term(dk);
    Sub     = term(ik);
    DT      = term(dkfac);  
    
    CHbl    = CT(keep,:);
    CHk     =  mean(CHbl,2);
    CH      =  term(CHk);
    
    M       = 1  + CH + A + S + G4 + DT + G4*DT + random(Sub) + I;
    slm = SurfStatLinModS(Ck,M,SW);
    
    
    %% Figure 2 
    for figure2  = 1
          
          contrast = ((((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) - ((G4.Control.* DT.T1)-(G4.Control.* DT.T0))))
          
          slm  	= SurfStatT(slm,contrast);
          statsG1   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
              [RPATH 'FIGURE2A.CT.G1.presence-control.T0T1.FAC.days.'],1);
          
          numclus = sum(statsG1.clus.P < 0.025);
          for i = 1 : numclus
              AAL2C2 = AAL(statsG1.clusid==i);
              unique(AAL2C2)
              vertexid = find(statsG1.t == max(statsG1.t(statsG1.clusid==i)));
              aalnum = AAL(:,vertexid);
              coordinate = SM.coord(:,vertexid)';
              maxtval = max(statsG1.t(statsG1.clusid==i));
              disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
              %      1     3     5     7     9    11    13    15    19    23    25    27    29    31    33
              % t-val = 5.0196 AAL: 3, coordinate: -15.1788      49.6029       36.329
              %      1     4     6     8    10    16    20    24    26    28    32    34    68
              % t-val = 5.072 AAL: 24, coordinate: 6.19002      27.0555      58.8839
          end
          
          numclus = sum(statsG1.inv_clus.P < 0.025);
          for i = 1 : numclus
              AAL2C2 = AAL(statsG1.inv_clusid==i);
              unique(AAL2C2)
              %         f= figure,
              %         hist(AAL(statsG1.inv_clusid==i),unique(AAL(statsG1.inv_clusid==i)));
              %         exportfigbo(f,[RPATH 'G2presenceN' num2str(i) '.png'],'png',10);
              vertexid = find(statsG1.inv_t == max(statsG1.inv_t(statsG1.inv_clusid==i)));
              aalnum = AAL(:,vertexid);
              coordinate = SM.coord(:,vertexid)';
              maxtval = max(statsG1.inv_t(statsG1.inv_clusid==i));
              disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
              %      1    51    53    55    57    61    63    65    81    85    89
              % t-val = 5.3919 AAL: 85, coordinate: -56.267     -45.6577     -12.1994
              %      1    40    43    44    46    48    50    52    54    68
              % t-val = 4.057 AAL: 44, coordinate: 7.81859     -96.7332      7.14824
              %      1    64    66    82    84    86    88    90
              % t-val = 4.3105 AAL: 86, coordinate: 57.2689     -2.11664     -31.0907
              %      2    12    14    18    58    64
              % t-val = 3.6607 AAL: 64, coordinate: 58.8944     -21.1868      27.0442
              %      1    43    45    47    49    51    53
              % t-val = 4.1444 AAL: 51, coordinate: -16.61529     -100.7017       4.65367
          end
          
          contrast = ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) - ((G4.Control.* DT.T1)-(G4.Control.* DT.T0))))
          slm  	= SurfStatT(slm,contrast);
          statsG2   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
              [RPATH 'FIGURE2A.CT.G2.presence-control.T0T1.FAC.days.'],1);
          
          numclus = sum(statsG2.clus.P < 0.025);
          for i = 1 : numclus
              AAL2C2 = AAL(statsG2.clusid==i);
              unique(AAL2C2)
              vertexid = find(statsG2.t == max(statsG2.t(statsG2.clusid==i)));
              aalnum = AAL(:,vertexid);
              coordinate = SM.coord(:,vertexid)';
              maxtval = max(statsG2.t(statsG2.clusid==i));
              disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
              %      1     4     6     8    10    14    16    20    24    26    28    32    34
              % t-val = 4.6543 AAL: 16, coordinate: 44.6538      48.6399     -7.64978
          end
          
          numclus = sum(statsG2.inv_clus.P < 0.025);
          for i = 1 : numclus
              AAL2C2 = AAL(statsG2.inv_clusid==i);
              unique(AAL2C2)
              %f= figure,
              %hist(AAL(statsG2.inv_clusid==i),unique(AAL(statsG2.inv_clusid==i)));
              %exportfigbo(f,[RPATH 'G2presenceN' num2str(i) '.png'],'png',10);
              vertexid = find(statsG2.inv_t == max(statsG2.inv_t(statsG2.inv_clusid==i)));
              aalnum = AAL(:,vertexid);
              coordinate = SM.coord(:,vertexid)';
              maxtval = max(statsG2.inv_t(statsG2.inv_clusid==i));
              disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
              %     40    44    46    48    68
              % t-val = 4.3292 AAL: 44, coordinate: 15.6613     -56.3666      8.03368
              %      1    18    30    64    80    82
              % t-val = 3.4476 AAL: 64, coordinate: 59.5053     -36.5102      23.0216
              %      7     9    13    15
              % t-val = 3.571 AAL: 13, coordinate: -40.1846      29.9264      13.0839
              
          end
           
          inter_tval = (statsG1.t) .* (statsG2.t)
          inter_len = (statsG1.p<p_thresh) .* (statsG2.p<p_thresh)
          f=figure;
          BoSurfStatViewData(double(inter_len>0),SInf,'');
          colormap([0 0 0; 1 0 0])
          exportfigbo(f,[RPATH 'FIGURE1.intersect_TC1-RCC_and_TC2-RCC_len.png'],'png',10);
          
          inter_neg_len = (statsG1.inv_p<p_thresh) .* (statsG2.inv_p<p_thresh)
          f=figure;
          BoSurfStatViewData(double(inter_neg_len>0),SM,'');
          colormap([0 0 0; 0 0 1])
          exportfigbo(f,[RPATH 'FIGURE1.intersect_TC1-RCC_and_TC2-RCC_invlen.png'],'png',10);
          
          % FIGURE1F -- OVERLAP
          inter_str = (statsG1.pval.C<pwe) .* (statsG2.pval.C<pwe)
          f=figure;
          BoSurfStatViewData(double(inter_str>0),SM,'');
          colormap([0 0 0; 1 0 0])
          exportfigbo(f,[RPATH 'FIGURE1.intersect_TC1-RCC_and_TC2-RCC_str.png'],'png',10);
          
          AAL2C2 = AAL(inter_str>0)
          unique(AAL2C2)
          hist(AAL(inter_str>0),unique(AAL(inter_str>0)))
          
          inter_neg_str = (statsG1.inv_pval.C<pwe) .* (statsG2.inv_pval.C<pwe)
          f=figure;
          BoSurfStatViewData(double(inter_neg_str>0),SM,'');
          colormap([0 0 0; 0 0 1])
          exportfigbo(f,[RPATH 'FIGURE1.intersect_TC1-RCC_and_TC2-RCC_invstr.png'],'png',10);
          
          AAL2C2 = AAL(inter_neg_str>0)
          unique(AAL2C2)
          hist(AAL(inter_neg_str>0),unique(AAL(inter_neg_str>0)))
         
          %% Functional mask
          load('attention_valcon.mat')
          
          f = figure,
          BoSurfStatViewData(o_attbin.*inter_len, SInf,'')
          exportfigbo(f,[RPATH 'ROI.att-bothlxF2.png'],'png',10);
          
          %% WHICH AAL does the overlap fall into?
          overlapF2 = o_attbin.*inter_len
          AAL2C2 = AAL(find(overlapF2>0));
          unique(AAL2C2)
          f= figure,
          hist(AAL(find(overlapF2>0)),unique(AAL(find(overlapF2>0))));
          % in order of size : 20, 24, 4, 8, 34, 51, 1
          
          f = figure,
          BoSurfStatViewData(o_attbin.*inter_neg_len, SInf,'')
          exportfigbo(f,[RPATH 'ROI.att-both.N.lxF1.png'],'png',10);
          
          overlapF2 = o_attbin.*inter_neg_len
          AAL2C2 = AAL(find(overlapF2>0));
          unique(AAL2C2)
          f= figure,
          hist(AAL(find(overlapF2>0)),unique(AAL(find(overlapF2>0))));
          % 64, 82, 86, 44
      end
    
    for figure3A = 1
          %% T0 -> T2
          contrast = ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T0)) - ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T0))
          slm  	= SurfStatT(slm,contrast);
          T2T0stats   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
              [RPATH 'FIGURE2A.CT.T2-T0.affect-perspective.FAC.days.'],1);
          T2T0stats   = BoSurfStatPlotAllStatsNew(slm, SInf, mask, p_thresh, 0.025, ...
              [RPATH 'FIGURE2A.CT.INF.T2-T0.affect-perspective.FAC.days.'],1);
          F2C_stats = stats;
          
          contrast = ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T0)) - ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T0))
          slm  	= SurfStatT(slm,contrast);
          
          %% test for complex model with random effects
          numclus = sum(T2T0stats.clus.P < 0.025);
          for i = 1 : numclus
              seedx       = mean(CT(keep,T2T0stats.clusid==i),2);
              seedxt      = term(seedx);
              
              M2       = 1 + CH + A + S + G4 + DT + G4*DT + Sub;
              slm1 = SurfStatLinMod(seedx,M2);
              slm1 = SurfStatT(slm1, contrast);
              disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
              % ROI: 1,t= 3.9804
              % ROI: 1,t= 4.2322
          end
          
          %% test for complex model with random effects
          numclus = sum(T2T0stats.inv_clus.P < 0.025);
          for i = 1 : numclus
              seedx       = mean(CT(keep,T2T0stats.inv_clusid==i),2);
              seedxt      = term(seedx);
              
              M2       = 1 + CH + A + S + G4 + DT + G4*DT + Sub;
              slm1 = SurfStatLinMod(seedx,M2);
              slm1 = SurfStatT(slm1, contrast);
              disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
              % ROI: 1,t= -4.7594
              % ROI: 1,t= -4.7681
          end
          
          
          %% Affect - perspective
          numclus = sum(T2T0stats.clus.P < 0.025);
          for i = 1 : numclus
              vertexid = find(T2T0stats.t == max(T2T0stats.t(T2T0stats.clusid==i)));
              aalnum = AAL(:,vertexid);
              coordinate = SM.coord(:,vertexid)';
              maxtval = max(T2T0stats.t(T2T0stats.clusid==i));
              AAL2C2 = AAL(T2T0stats.clusid==i);
              unique(AAL2C2)
              disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
              %      1     3     7    13    19    23    31    33
              % t-val = 3.8197 AAL: 23, coordinate: -12.0624      36.2348      22.4073
              %      1    19    33    67    69
              % t-val = 3.4648 AAL: 33, coordinate: -11.2437     -3.47004       43.516
          end
          
          %% Perspective - affect
          numclus = sum(T2T0stats.inv_clus.P < 0.025);
          for i = 1 : numclus
              vertexid = find(T2T0stats.inv_t == max(T2T0stats.inv_t(T2T0stats.inv_clusid==i)));
              aalnum = AAL(:,vertexid);
              coordinate = SM.coord(:,vertexid)';
              maxtval = max(T2T0stats.inv_t(T2T0stats.inv_clusid==i));
              AAL2C2 = AAL(T2T0stats.inv_clusid==i);
              unique(AAL2C2)
              disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
              % 1    17    57    61    63    65    81    85    89
              % t-val = 4.2654 AAL: 81, coordinate: -64.3335     -39.6449        13.92
              % 1    82    84    86    88    90
              % t-val = 4.9318 AAL: 86, coordinate: 56.1768     -16.6383     -15.5999
          end
          
          %% Functional masks
          load('empatom_tom.mat')
          load('empatom_com.mat')
          
          %% Figure 3 A
          F2tombin = T2T0stats.inv_p<0.025
          f = figure,
          BoSurfStatViewData(o_etombin.*F2tombin, SInf,'')
          exportfigbo(f,[RPATH 'ROI.ToMFA.png'],'png',10);
          %% WHICH AAL does the overlap fall into?
          overlapF2 =o_etombin.*F2tombin
          AAL2C2 = AAL(find(overlapF2>0));
          unique(AAL2C2)
          f= figure,
          hist(AAL(find(overlapF2>0)),unique(AAL(find(overlapF2>0))));
          % 81, 63, 85, 82, 86, 61, 88, 65, 13, 17
          
          F2affbin = T2T0stats.p<0.025
          f = figure,
          BoSurfStatViewData(o_eempbin.*F2affbin, SInf,'')
          exportfigbo(f,[RPATH 'ROI.AffFA.png'],'png',10);
          
          overlapF2 = o_eempbin.*F2affbin
          AAL2C2 = AAL(find(overlapF2>0));
          unique(AAL2C2)
          f= figure,
          hist(AAL(find(overlapF2>0)),unique(AAL(find(overlapF2>0))));
          % 23, 24, 68, 67, 8, 13, 19, 3, 4, 7, 12, 14, 15, 20, 29, 30
          
      end
    
    for figure3B = 1
          % (Perspective - Affect ) - Presence
          contrast_per = (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) + ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)))
          contrast_aff = (((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) + ((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)))
          contrast_pre = (((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_2.* DT.T0)))
          contrast_per_vs_aff_vs_pre = (contrast_per-contrast_aff)-contrast_pre;
          
          slm  	= SurfStatT(slm,contrast_per_vs_aff_vs_pre);
          perspectiveF2.stats   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
              [RPATH 'FIGURE2A.CT.perspective-affect.FAC.days.M3.'],1);
          perspectiveF2.stats   = BoSurfStatPlotAllStats(slm, SInf, mask, p_thresh, 0.025, ...
              [RPATH 'FIGURE2A.CT.inf.perspective-affect.FAC.days.'],1);
          
          %% test for complex model with random effects // FIXED effects
          numclus = sum(perspectiveF2.stats.clus.P < 0.025);
          for i = 1 : numclus
              seedx       = mean(CT(keep,perspectiveF2.stats.clusid==i),2);
              seedxt      = term(seedx);
              
              M2       = 1 + CH + A + S + G4 + DT + G4*DT + Sub;
              slm1 = SurfStatLinModS(seedx,M2);
              slm1 = SurfStatT(slm1, contrast_per_vs_aff_vs_pre);
              disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
              % ROI: 1,t= 4.0308
              % ROI: 1,t= 4.4525
              % ROI: 1,t= 3.7377
              % ROI: 1,t= 3.5254
              %         M3       = CH + A + S + (1 + G4)*(1 + DT)*(1 + random(Sub)) + I;
              %         slm1 = SurfStatLinMod(seedx,M3);
              %         slm1 = SurfStatT(slm1, contrast);
              %         disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
          end
          
          numclus = sum(perspectiveF2.stats.clus.P < 0.025);
          for i = 1 : numclus
              vertexid = find(perspectiveF2.stats.t == max(perspectiveF2.stats.t(perspectiveF2.stats.clusid==i)));
              aalnum = AAL(:,vertexid);
              coordinate = SM.coord(:,vertexid)';
              maxtval = max(perspectiveF2.stats.t(perspectiveF2.stats.clusid==i));
              AAL2C2 = AAL(perspectiveF2.stats.clusid==i);
              unique(AAL2C2)
              disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
              %      1    81    85    89
              % t-val = 3.6375 AAL: 85, coordinate: -51.2848     -43.3906     0.232559
              %      1    82    86    90
              % t-val = 4.0353 AAL: 86, coordinate: 65.302     -42.8761     -8.28508
              %      1    43    44    46    48    50
              % t-val = 3.0493 AAL: 48, coordinate: 10.4429     -84.4072     -13.6755
              %      1    43    45    47    49    51
              % t-val = 3.0419 AAL: 43, coordinate: -10.8362     -86.9299      7.04825
          end
          
          
          %% (Affect - Perspective) - Presence
          contrast_aff_vs_per_vs_pre = (contrast_aff-contrast_per)-contrast_pre;
          slm  	= SurfStatT(slm,contrast_aff_vs_per_vs_pre);
          affectiveF2.stats   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
              [RPATH 'FIGURE2A.CT.affect-perspective.FAC.M.days.'],1);
          affectiveF2.stats   = BoSurfStatPlotAllStats(slm, SInf, mask, p_thresh, 0.025, ...
              [RPATH 'FIGURE2A.CT.INF.affect-perspective.b.FAC.days.'],1);
          
          numclus = sum(affectiveF2.stats.clus.P < 0.025);
          for i = 1 : numclus
              seedx       = mean(CT(keep,affectiveF2.stats.clusid==i),2);
              seedxt      = term(seedx);
              
              M2       = 1 + CH + A + S + G4 + DT + G4*DT + Sub;
              slm1 = SurfStatLinMod(seedx,M2);
              slm1 = SurfStatT(slm1, contrast_aff_vs_per_vs_pre);
              disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
              % ROI: 1,t= 5.6359
              % ROI: 1,t= 4.1286
              % ROI: 1,t= 3.7002
              % ROI: 1,t= 3.1319
          end
          
          %% AAL
          numclus = sum(affectiveF2.stats.clus.P < 0.025);
          for i = 1 : numclus
              vertexid = find(affectiveF2.stats.t == max(affectiveF2.stats.t(affectiveF2.stats.clusid==i)));
              aalnum = AAL(:,vertexid);
              coordinate = SM.coord(:,vertexid)';
              maxtval = max(affectiveF2.stats.t(affectiveF2.stats.clusid==i));
              AAL2C2 = AAL(affectiveF2.stats.clusid==i);
              unique(AAL2C2)
              disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
              %      1     2     8    12    14    16    18    30    58    64
              %      t-val = 4.1065 AAL: 14, coordinate: 53.3627      27.3439       1.9629
              %      1    33    35    67    69
              %      t-val = 3.4612 AAL: 33, coordinate: -7.97467     -40.1174      38.1422
              %      1    61    63    65    85
              %      t-val = 3.3241 AAL: 61, coordinate: -51.5047     -50.0944      36.4772
              %      1    51    53    55    85    89
              %      t-val = 2.9988 AAL: 85, coordinate: -63.1341     -38.1133      -14.029
          end
          
          % Functional masks
          
          %% Figure 3B
          F2Btombin = perspectiveF2.stats.puncorr<0.025
          f = figure,
          BoSurfStatViewData(o_etombin.*F2Btombin, SInf,'')
          exportfigbo(f,[RPATH 'ROI.ToMFB.png'],'png',10);
          
          overlapF2 = o_etombin.*F2Btombin
          AAL2C2 = AAL(find(overlapF2>0));
          unique(AAL2C2)
          f= figure,
          hist(AAL(find(overlapF2>0)),unique(AAL(find(overlapF2>0))));
          % 81, 13, 85, 82, 86
          
          F2Baffbin = affectiveF2.stats.puncorr<0.025
          f = figure,
          BoSurfStatViewData(o_eempbin.*F2Baffbin, SInf,'')
          exportfigbo(f,[RPATH 'ROI.AffFB.png'],'png',10);
          
          overlapF2 = o_eempbin.*F2Baffbin
          AAL2C2 = AAL(find(overlapF2>0));
          unique(AAL2C2)
          f= figure,
          hist(AAL(find(overlapF2>0)),unique(AAL(find(overlapF2>0))));
          % 14, 61, 65, 67, 8, 12, 33, 35, 63, 16, 2, 64
          
          F2Btombin = perspectiveF2.stats.puncorr<0.025
          f = figure,
          BoSurfStatViewData(F2Btombin, SInf,'')
          exportfigbo(f,[RPATH 'ROI.PFB.png'],'png',10);
          
          F2Baffbin = affectiveF2.stats.puncorr<0.025
          f = figure,
          BoSurfStatViewData(F2Baffbin, SInf,'')
          exportfigbo(f,[RPATH 'ROI.AFB.png'],'png',10);
      end
      
    %% supplements
    
    for subfigure1 = 1
        %% presence versus perspective
        presence_contrast    = (((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)))
        perspective_contrast = (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) + ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)))
        between_affect       = ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1))
        contrast = (presence_contrast - (perspective_contrast - between_affect))
        slm  	   = SurfStatT(slm,contrast);
        statsprp   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.presence-perspective.v3.FAC.days.'],1);
        
        contrast = ((perspective_contrast - between_affect) - presence_contrast)
        slm  	   = SurfStatT(slm,contrast);
        statsppr   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.perspective-presence.v3.FAC.days.'],1);
        
        numclus = sum(statsprp.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(statsprp.clusid==i);
            unique(AAL2C2)

            vertexid = find(statsprp.t == max(statsprp.t(statsprp.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(statsprp.t(statsprp.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %    1     2     4     6     8    10    12    14    16    18    20    24    30
            % t-val = 4.606 AAL: 8, coordinate: 22.6024      25.5355      42.3591
        end
        
        
        numclus = sum(statsppr.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(statsppr.clusid==i);
            unique(AAL2C2)

            vertexid = find(statsppr.t == max(statsppr.t(statsppr.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(statsppr.t(statsppr.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            
            %    1     3     5     7     9    11    13    15    23    25
            % t-val = 4.6562 AAL: 13, coordinate: -43.8494      27.9434      12.6475
            %      1    43    45    47    49
            % t-val = 3.3158 AAL: 47, coordinate: -8.57926     -85.0927     -11.3089
        end
        
        
        %% presence versus affect
        contrast   = ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)))) - ...
            ((((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)) + ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1))) - ...
            (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1))))
        
        presence_contrast = (((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)))
        affect_contrast = (((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) + ((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)))
        between_perspective = ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1))
        contrast = (presence_contrast - (affect_contrast - between_perspective))
        slm  	   = SurfStatT(slm,contrast);
        statspra   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.presence-affect.v3.FAC.days.'],1);
        
        contrast   = ((((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)) + ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)))) - ...
            ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))) - ...
            (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1))))
        
        contrast = ((affect_contrast - between_perspective) - presence_contrast)
        slm  	   = SurfStatT(slm,contrast);
        statsapr   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.affect-presence.v3.FAC.days.'],1);
        
        
        
        numclus = sum(statspra.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(statspra.clusid==i);
            unique(AAL2C2)
            %         f= figure,
            %         hist(AAL(statspra.clusid==i),unique(AAL(statspra.clusid==i)));
            %         exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
            vertexid = find(statspra.t == max(statspra.t(statspra.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(statspra.t(statspra.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %      1     3     5     7     9    13    15    23    25    27
            % t-val = 3.6021 AAL: 7, coordinate: -31.7867      51.7903      2.16786
            %      1     6    10    16    26    28
            % t-val = 3.7732 AAL: 28, coordinate: 3.62128      42.8642     -25.1934
        end
        
        
        
        numclus = sum(statsapr.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(statsapr.clusid==i);
            unique(AAL2C2)
            f= figure,
            hist(AAL(statsapr.clusid==i),unique(AAL(statsapr.clusid==i)));
            exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
            vertexid = find(statsapr.t == max(statsapr.t(statsapr.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(statsapr.t(statsapr.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %    1     2     8    12    14    18    30    58    64
            % t-val = 3.731 AAL: 2, coordinate: 47.1103      2.52934      46.6577
            %      1    19    33    35    67    69
            % t-val = 3.8004 AAL: 33, coordinate: -7.97467     -40.1174      38.1422
        end
        
        
    end   % modules versus presence
    
    for subfigure2 = 2
        %% affect - controls
        contrastaff1   = (((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)) - ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) ...
            - ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)));
        contrastaff2   = (((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)));
        
        contrastcontrols = (((G4.Control.* DT.T3)-(G4.Control.* DT.T2)) + ((G4.Control.* DT.T2)-(G4.Control.* DT.T1)));
        contrast = (contrastaff1 + contrastaff2) - contrastcontrols;
        slm  	= SurfStatT(slm,contrast);
        statsa   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.affect-C.FAC.days.'],1);
        
        numclus = sum(statsa.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(statsa.clusid==i);
            unique(AAL2C2)
            
            vertexid = find(statsa.t == max(statsa.t(statsa.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(statsa.t(statsa.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %         1     2     8    12    14    18    30    58    64
            %         t-val = 4.1862 AAL: 2, coordinate: 49.5661     -1.86254      49.9535
            %         1    19    33    35    67    69
            %         t-val = 3.8925 AAL: 33, coordinate: -7.97467     -40.1174      38.1422
        end
        
        
        %% perspective - controls
        contrastpers1   = (((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)) - ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) ...
            - ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)));
        contrastpers2   = (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) - ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)));
        
        contrastcontrols = (((G4.Control.* DT.T3)-(G4.Control.* DT.T2)) + ((G4.Control.* DT.T2)-(G4.Control.* DT.T1)));
        contrast = (contrastpers1 + contrastpers2) - contrastcontrols;
        slm  	= SurfStatT(slm,contrast);
        statsp   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.perspective-Controls.FAC.days.'],1);
        
        numclus = sum(statsp.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(statsp.clusid==i);
            unique(AAL2C2)
            
            vertexid = find(statsp.t == max(statsp.t(statsp.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(statsp.t(statsp.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %         1     3     5     7     9    11    13    15    23    25    31
            %         t-val = 4.6463 AAL: 13, coordinate: -41.1369      27.9076      14.5408
        end
        
        contrast = ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) - ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0))
        slm  	= SurfStatT(slm,contrast);
        statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.perspectiveG2.FAC.days.'],1);
        
        contrast = ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)) - ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ...
            ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))
        slm  	= SurfStatT(slm,contrast);
        statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.perspectiveG1.FAC.days.'],1);
        
    end    % modules versus controls
    
    for subfigure3 = 3       % use for overlap of CTX and DBM
        
        contrast = ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))))- ...
            ((G4.Control.* DT.T1)-(G4.Control.* DT.T0))
        slm  	= SurfStatT(slm,contrast);
        statsG1G2   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.G1G2.presence-control.T0T1.FAC.days.'],1);
    end

    %% not used now but maybe keep it in for now: T1-T2 and T2-T3
    for notused = 1
        %% S1: T1_T2:
        contrast_T1_T2    = (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) + ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T1)))
        contrast_presence = (((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)))
        contrast =  contrast_T1_T2 - contrast_presence;
        
        
        slm  	= SurfStatT(slm,contrast);
        statsG2G1_T21   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.G2-G1.T1T2.FAC.days.'],1);
        
        
        numclus = sum(statsG2G1_T21.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(statsG2G1_T21.clusid==i);
            unique(AAL2C2)
            %f= figure,
            %hist(AAL(statsG2G1_T21.clusid==i),unique(AAL(statsG2G1_T21.clusid==i)));
            %exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
            vertexid = find(statsG2G1_T21.t == max(statsG2G1_T21.t(statsG2G1_T21.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(statsG2G1_T21.t(statsG2G1_T21.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %      1     3     5     7     9    11    13    15    19    23    25    27    29    31    33
            % t-val = 6.1317 AAL: 7, coordinate: -36.6672      54.3752      2.32347
        end
        
        
        contrast_T1_T2_12    = (((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)))
        contrast_presence    = (((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) - ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)))
        contrast =  contrast_T1_T2_12 - contrast_presence;
        
        slm  	= SurfStatT(slm,contrast);
        statsG1G2_T21   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.G1-G2.T1T2.FAC.days.'],1);
        
        
        numclus = sum(statsG1G2_T21.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(statsG1G2_T21.clusid==i);
            unique(AAL2C2)
            %         f= figure,
            %         hist(AAL(statsG1G2_T21.clusid==i),unique(AAL(statsG1G2_T21.clusid==i)));
            %         exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
            vertexid = find(statsG1G2_T21.t == max(statsG1G2_T21.t(statsG1G2_T21.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(statsG1G2_T21.t(statsG1G2_T21.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            
            % 1    51    53    55    61    63    65    81    85    89
            % t-val = 4.5286 AAL: 61, coordinate: -57.5676     -44.5012      39.0814
            % 1     2     4     8    10    12    14    16    18    30    58    64
            % t-val = 5.5108 AAL: 14, coordinate: 53.3548      31.2898      3.31863
            % 82    84    86    88    90
            % t-val = 4.0061 AAL: 86, coordinate: 60.0739     0.359326     -26.7405
        end
        
        
        %% T3_T2:
        
        contrast = ((((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)) - ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)))) - ...
            (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) - ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1))) - ...
            (((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) - ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)))
        
        contrast_T2_T3_21    = (((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)) - ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)))
        contrast_T1_T2_21    = (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) - ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)))
        contrast_presence    = (((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) - ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)))
        contrast =  contrast_T2_T3_21 - contrast_T1_T2_21 - contrast_presence;
        
        slm  	= SurfStatT(slm,contrast);
        statsG2G1_T23   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.G2-G1.T3T2+.FAC.days.'],1);
        
        
        
        numclus = sum(statsG2G1_T23.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(statsG2G1_T23.clusid==i);
            unique(AAL2C2)
            %         f= figure,
            %         hist(AAL(statsG2G1_T23.clusid==i),unique(AAL(statsG2G1_T23.clusid==i)));
            %         exportfigbo(f,[RPATH 'G2_G1_T2T3' num2str(i) '.png'],'png',10);
            vertexid = find(statsG2G1_T23.t == max(statsG2G1_T23.t(statsG2G1_T23.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(statsG2G1_T23.t(statsG2G1_T23.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            
            %      1    19    33    67    69
            % t-val = 3.655 AAL: 33, coordinate: -16.5865     -32.0004      43.0852
            %      2     4     8    24
            % t-val = 3.2665 AAL: 8, coordinate: 41.9548      12.7887      50.0459
            %     39    47    55
            % t-val = 4.1168 AAL: 39, coordinate: -17.0702     -37.5658     -13.4916
        end
        
        
        contrast = ((((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)) - ((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)))) - ...
            (((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1))) - ...
            (((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) - ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)))
        
        contrast_T2_T3_12    = (((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)) - ((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)))
        contrast_T1_T2_12    = (((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)))
        contrast_presence    = (((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) - ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)))
        contrast =  contrast_T2_T3_12 - contrast_T1_T2_12 - contrast_presence;
        
        slm  	= SurfStatT(slm,contrast);
        statsG1G2_T23  = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.CT.G1-G2.T3T2+.FAC.days.'],1);
        
        
        numclus = sum(statsG1G2_T23.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(statsG1G2_T23.clusid==i);
            unique(AAL2C2)
            %         f= figure,
            %         hist(AAL(statsG1G2_T23.clusid==i),unique(AAL(statsG1G2_T23.clusid==i)));
            %         exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
            vertexid = find(statsG1G2_T23.t == max(statsG1G2_T23.t(statsG1G2_T23.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(statsG1G2_T23.t(statsG1G2_T23.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            
            %      1    17    57    61    63    81    85
            % t-val = 4.3501 AAL: 85, coordinate: -52.8215     -48.3212      3.25627
            %      1    82    84    86    88    90
            % t-val = 4.5855 AAL: 86, coordinate: 56.1768     -16.6383     -15.5999
        end
        
    end
     
end

%% DBM on SURFACE CHANGE | SUPPLEMENTs
for SUB = 1
    keep    = mintersect(find(isthere_dbm), ...
        find(tp_t3_there==1), ...
        find(tp_t1_there==1), ...
        find(tp_t2_there==1 ), ...
        find(tp_t0_there==1 ), ...
        union(find(tpnum==0), union(find(tpnum==1),union(find(tpnum==3), find(tpnum==2 )))));
    
    mask    = mean(CT,1)>0.4;
    Ck      = zdbm_map(keep,:);
    ik      = id(keep,:);
    ink     = idnum(keep);
    tk      = tpnum(keep,:);
    dk      = Days_last(keep);
    dkfac   = tp(keep);
    gk      = group(keep);
    g4k     = group4(keep);
    ak      = age(keep);
    sk      = sex(keep);
    
    A       = term(ak);
    S       = term(sk);
    G       = term(gk);
    G4      = term(g4k);
    D       = term(dk);
    Sub     = term(ik);
    DT      = term(dkfac);
    
    CHk =  mean(Ck,2);
    CH  =  term(CHk);

    M       = 1 + CH + A + S + G4 + DT + G4*DT + Sub;
    slm     = SurfStatLinModS(Ck,M,SW);

    for subfigure = 1
        % presence 
        contrast = (((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))) - ...
            (((G4.Control.* DT.T1)-(G4.Control.* DT.T0)))
        slm  	= SurfStatT(slm,contrast);
        stats   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE1A.dbm.group1group2-RCC.FAC.days.'],1);
        dbm_presece = stats;
        
        numclus = sum(stats.clus.P < 0.025);
        for i = 1 : numclus
            vertexid = find(stats.t == max(stats.t(stats.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(stats.t(stats.clusid==i));
            AAL2C2 = AAL(stats.clusid==i);
            unique(AAL2C2)
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %         4     6     8    10    14    16    24    26
            %         t-val = 5.056 AAL: 26, coordinate: 9.58192      66.1969     -7.40471
            %         45    49    51    57    59    61    65    67    69
            %         t-val = 4.153 AAL: 59, coordinate: -29.7819     -57.4046       59.997
            %         3     5     7     9    13    15    23    25
            %         t-val = 5.2948 AAL: 9, coordinate: -39.4862      47.7216     -5.30296
            %         50    60    68
            %         t-val = 4.3725 AAL: 60, coordinate: 18.0427     -69.8457      48.3119
        end
        
        numclus = sum(stats.inv_clus.P < 0.025);
        for i = 1 : numclus
            vertexid = find(stats.inv_t == max(stats.inv_t(stats.inv_clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(stats.inv_t(stats.inv_clusid==i));
            AAL2C2 = AAL(stats.inv_clusid==i);
            unique(AAL2C2)
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            
            %      1    43    45    47    51    53    55    67    85    89
            % t-val = 6.5234 AAL: 47, coordinate: -6.14305     -67.2225    -0.316041
            %      1    43    44    46    48    54    56    68
            % t-val = 6.2539 AAL: 48, coordinate: 12.229     -71.7176     -7.65261
            %     82    84    86    88    90
            % t-val = 4.3318 AAL: 90, coordinate: 55.9764     -13.7031     -26.6315
            %      1    57    61    63
            % t-val = 4.204 AAL: 1, coordinate: -59.3612     0.400031      30.3988
            %     64    66    82
            % t-val = 3.8853 AAL: 82, coordinate: 59.9256     -13.6744       2.8411
            %     81    83    85
            % t-val = 4.0227 AAL: 85, coordinate: -49.7752    -0.618724     -23.5408
            %      7    13
            % t-val = 3.7831 AAL: 13, coordinate: -55.2127      21.9125      17.2652
            %     66    82    86    90
            % t-val = 3.9725 AAL: 86, coordinate: 52.3621     -58.1651      10.0281
        end
        
          %% DBM and CT
          CTF1 = statsG1G2.p<0.025
          PCF1 = dbm_presece.p<0.025
          f = figure,
          BoSurfStatViewData(CTF1.*PCF1, SInf,'')
          exportfigbo(f,[RPATH 'presence.overlapCTDBM.png'],'png',10);
          CTF1 = statsG1G2.inv_p<0.025
          PCF1 = dbm_presece.inv_p<0.025
          f = figure,
          BoSurfStatViewData(CTF1.*PCF1, SInf,'')
          exportfigbo(f,[RPATH 'presence.overlapCTDBM.N.png'],'png',10);

        
        % for B T2-T0
        contrast = (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T0)) - ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T0)))
        slm  	= SurfStatT(slm,contrast);
        stats   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE1A.dbm.TC2TC1.T0T2.FAC.days.'],1);
        dbm_T0T2_stats21 = stats;
        
        numclus = sum(stats.clus.P < 0.025);
        for i = 1 : numclus
            vertexid = find(stats.t == max(stats.t(stats.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(stats.t(stats.clusid==i));
            AAL2C2 = AAL(stats.clusid==i);
            unique(AAL2C2)
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %      1    81    85
            % t-val = 4.8892 AAL: 1, coordinate: -48.2206     -32.9208     -7.11724
            %      1    51    61    63    65    81    85
            % t-val = 4.6725 AAL: 85, coordinate: -53.6574     -55.8975      21.4058
            %     64    82    84    86    88
            % t-val = 3.4369 AAL: 82, coordinate: 61.5434     -3.36827    -0.291067
            %     11    13    17    29    57
            % t-val = 3.4516 AAL: 11, coordinate: -44.3014      12.7902      1.40068
        end
        
        numclus = sum(stats.inv_clus.P < 0.025);
        for i = 1 : numclus
            vertexid = find(stats.inv_t == max(stats.inv_t(stats.inv_clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(stats.inv_t(stats.inv_clusid==i));
            AAL2C2 = AAL(stats.inv_clusid==i);
            unique(AAL2C2)
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %      1     2     4     6     8    10    12    14    18    20    24    26    28    32    34    58    60    62    68    70
            % t-val = 5.8279 AAL: 26, coordinate: 8.37617      56.2547     -7.43806
            %      1     3     5     7     9    15    19    23    25    31    33    57    59    67    69
            % t-val = 4.2862 AAL: 3, coordinate: -24.691      59.5709      14.7184
        end
        
          %% T2-T0
          CTF1 = T2T0stats.p<0.025
          PCF1 = dbm_T0T2_stats21.inv_p<0.025
          f = figure,
          BoSurfStatViewData(PCF1.*CTF1, SInf,'')
          exportfigbo(f,[RPATH 'affectT0T2.overlapCTDBM.png'],'png',10);
          CTF1 = T2T0stats.inv_p<0.025
          PCF1 = dbm_T0T2_stats21.p<0.025
          f = figure,
          BoSurfStatViewData(PCF1.*CTF1, SInf,'')
          exportfigbo(f,[RPATH 'perspectiveT0T2.overlapCTDBM.png'],'png',10);
  
        
        % for C (perspective versus affective)
        contrast_per = (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) + ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)))
        contrast_aff = (((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) + ((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)))
        contrast_pre = (((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_2.* DT.T0)))
        contrast_per_vs_aff_vs_pre = (contrast_per-contrast_aff)-contrast_pre;
        contrast_aff_vs_per_vs_pre = (contrast_aff-contrast_per)-contrast_pre;
        
        %
        slm  	= SurfStatT(slm,contrast_aff_vs_per_vs_pre);
        stats   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.dbm.affect-perspective.FAC.days.'],1);
        dbm_ap_stats = stats;
        
        numclus = sum(stats.clus.P < 0.025);
        for i = 1 : numclus
            vertexid = find(stats.t == max(stats.t(stats.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(stats.t(stats.clusid==i));
            AAL2C2 = AAL(stats.clusid==i);
            unique(AAL2C2)
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            
            %         61    63    65
            %         t-val = 4.3277 AAL: 61, coordinate: -51.3293      -50.365       39.036
            %         57    61    63
            %         t-val = 3.7878 AAL: 57, coordinate: -59.9692     -8.42342      29.6387
        end
        
        slm  	= SurfStatT(slm,contrast_per_vs_aff_vs_pre);
        stats   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
            [RPATH 'FIGURE2A.dbm.perspective-affect.FAC.days.'],1);
        dbm_pa_stats = stats;
        
        numclus = sum(stats.clus.P < 0.025);
        for i = 1 : numclus
            clusadapt = ((stats.clusid==i));
            tvaluesclus = stats.t(clusadapt);
            maskt   = tvaluesclus<10; % made this for a weird outlier
            vertexid = find(stats.t == max(tvaluesclus(maskt)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max((tvaluesclus(maskt)));
            AAL2C2 = AAL(clusadapt);
            unique(AAL2C2)
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %      1    43    47    51    53    55
            % t-val = 5.3546 AAL: 47, coordinate: -24.4063      -73.736      -5.1867
            %      1    81    85
            % t-val = 3.6909 AAL: 85, coordinate: -51.5811     -10.4612     -15.4598
            %      1    40    56    84    88    90
            % t-val = 3.84 AAL: 84, coordinate: 40.80      7.11822     -24.9645
            %     82    84    86    88
            % t-val = 3.4898 AAL: 88, coordinate: 45.7671      9.66348     -30.1877
        end
        
        %% DIFF
        CTF1 = affectiveF2.stats.puncorr<0.025
        PCF1 = dbm_ap_stats.puncorr<0.025
        f = figure,
        BoSurfStatViewData(PCF1.*CTF1, SInf,'')
        exportfigbo(f,[RPATH 'affect.overlapCTDBM.png'],'png',10);
        CTF1 = perspectiveF2.stats.puncorr<0.025
        PCF1 = dbm_pa_stats.puncorr<0.025
        f = figure,
        BoSurfStatViewData(PCF1.*CTF1, SInf,'')
        exportfigbo(f,[RPATH 'perspective.overlapCTDBM.png'],'png',10);
        
    end

end

%% behavior
for whotokeep = 1
    keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
    keep1   = union(keep1,    find(strcmp(groupN,'Presence')));
    %keep1   = union(keep1,    find(strcmp(groupN,'Control')));
    
    keep4   = find(~(strcmp(group4,'Group3')));
    keep2   = intersect(find(X(:,1) >-666), find(sum(X,2)~=0));
    size(mintersect(keep1,keep2,keep4)) % N = 537 // 339 without controls
end




%% Compassion x AFFECT   
for behavioralaffect = 1
    for figurec = 1
        keep1   = find((strcmp(groupN,'Affect')))
        keep2   = intersect(find(X(:,1) >-666), find(sum(X,2)~=0));
        keep3   = find(abs(mean(X,2))<0.1);
        keep4   = find(~(strcmp(group4,'Group3')));
        size(mintersect(keep1,keep2,keep3, keep4)) % N = 126 (127)
        mask    = mean(CT,1)>0.9;
        keep7   = find(compassion>0);
        keep7b  = find(abs(Comp)<666);
        keep    = mintersect(keep1,keep2,keep4,keep7,keep7b, keep3);

        size(keep)
        size(mintersect(keep, find(strcmp(groupN,'Affect')))) % N = 121
        
        Ck      = X(keep,:);
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = sex(keep);
        Compk   = (Comp(keep));
        dk      = Days_last(keep);
        combf   = (chcbl(keep));
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ik);
        CMP     = term(Compk);
        T       = term(tk);
        Tn      = term(tnk);
        D       = term(dk);
        
        CHk   = mean(Ck,2);
        CH    = term(CHk);
        
        Cbf   = term(combf);
        
        M       = 1 + CH + A + S + CMP + random(Sub) + I;
        slm     = SurfStatLinModS(Ck,M,SW);
        slm  	= SurfStatT(slm,(Compk));
        
        f=figure,
        BoSurfStatViewData(slm.t,SInf,'')
        SurfStatColLim([-3 3])
        colormap([0 0 0; mycol.blackblue])
        exportfigbo(f,[RPATH 'FIG3C.CMP_blank.png'],'png',10)
        cluscomp   = BoSurfStatPlotAllStatsNew(slm,SM, mask, 0.025, 0.025, [RPATH '/FIG3C.affect.comp.WB' ],1);
        cluscomp   = BoSurfStatPlotAllStats(slm,SInf, mask, 0.025, 0.025, [RPATH '/FIG3C.affect.INF.comp.WB' ],1);
        
        
        %% test for complex model with random effects
        numclus = sum(cluscomp.clus.P < 0.025);
        for i = 1 : numclus
            seedx       = mean(X(keep,cluscomp.clusid==i),2);
            seedxt      = term(seedx);
            
            M2       = 1 + CH + A + S + GN + CMP + GN*CMP + Sub;
            slm1 = SurfStatLinMod(seedx,M2);
            slm1 = SurfStatT(slm1, (GN.Affect.*Compk));
            disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
            
            M3       = CH + A + S + (1 + CMP)*(1 + GN)*(1 + random(Sub)) + I;
            slm1 = SurfStatLinModS(seedx,M3);
            slm1 = SurfStatT(slm1, (GN.Affect.*Compk));
            disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
            
            % ROI: 1,t= 3.05
            % ROI: 1,t= 2.04
        end
        
        
        
        %which AAL?
        numclus = sum(cluscomp.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(cluscomp.clusid==i);
            unique(AAL2C2)
            f= figure,
            hist(AAL(cluscomp.clusid==i),unique(AAL(cluscomp.clusid==i)));
            exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
            vertexid = find(cluscomp.t == max(cluscomp.t(cluscomp.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(cluscomp.t(cluscomp.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %NEW
            %     1    12    14    16    30
            % t-val = 3.9006 AAL: 14, coordinate: 32.1758      22.3563      9.10549
            
        end
        
        numclus = sum(cluscomp.clus.P<0.025)
        for j = [1 : numclus]
            target = mean(Ck(:,find(cluscomp.clusid==j)),2);
            M    = 1 + CH ;
            slm = SurfStatLinMod(target,M);
            targetresid(keep,:) = target - slm.X*slm.coef;
            keepc1  = find(abs(Comp)<666);
            keepp = find(strcmp(groupN,'Affect'));
            keepc = mintersect(keep,keepc1,keepp);
            f= figure,
            hold on
            scatter(Comp(keepc),targetresid(keepc),60,'Marker','o', 'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0 0 0]), lsline
            %axis([-1.5 2.5 -0.65 0.4])
            exportfigbo(f,[RPATH 'FIG3C.Comp.affect.scatter' num2str(j) '.png'],'png',10)
            [r p] = corr(targetresid(keepc),Comp(keepc))
            %r =0.23
        end
        
        %% separate groups
        numclus = sum(cluscomp.clus.P<0.025)
        for j = [1 : numclus]
            target = mean(Ck(:,find(cluscomp.clusid==j)),2);
            M    = 1 + CH;
            slm = SurfStatLinMod(target,M);
            targetresid(keep,:) = target - slm.X*slm.coef;
            keepc1  = find(abs(Comp)<666);
            keepp = find(strcmp(groupN,'Affect'));
            keepp2 = find(strcmp(group,'Group_1'));
            keepc = mintersect(keep,keepc1,keepp,keepp2);
            f= figure,
            hold on
            scatter(Comp(keepc),targetresid(keepc),60,'Marker','o', 'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0 0 0]), lsline
            axis([-1.5 2.5 -0.2 0.2])
            exportfigbo(f,[RPATH 'FIG3C.Comp.affect.scatter.G1.' num2str(j) '.png'],'png',10)
            [r p] = corr(targetresid(keepc),Comp(keepc))
            %r =0.31, p<0.02
        end
        
        numclus = sum(cluscomp.clus.P<0.025)
        for j = [1 : numclus]
            target = mean(Ck(:,find(cluscomp.clusid==j)),2);
            M    = 1 + CH;
            slm = SurfStatLinMod(target,M);
            targetresid(keep,:) = target - slm.X*slm.coef;
            keepc1  = find(abs(Comp)<666);
            keepp = find(strcmp(groupN,'Affect'));
            keepp2 = find(strcmp(group,'Group_2'));
            keepc = mintersect(keep,keepc1,keepp,keepp2);
            f= figure,
            hold on
            scatter(Comp(keepc),targetresid(keepc),60,'Marker','o', 'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0 0 0]), lsline
            axis([-1.5 2.5 -0.2 0.2])
            exportfigbo(f,[RPATH 'FIG3C.Comp.affect.scatter.G2.' num2str(j) '.png'],'png',10)
            [r p] = corr(targetresid(keepc),Comp(keepc))
            %r =0.21
        end
        
            %% compassion x affect
    
    comF3 = cluscomp.puncorr<0.025
    f = figure,
    BoSurfStatViewData(o_eempbin.*comF3, SInf,'')
    exportfigbo(f,[RPATH 'ROI.com*aff.png'],'png',10);
    
    overlapF2 = o_eempbin.*comF3
    AAL2C2 = AAL(find(overlapF2>0));
    unique(AAL2C2)
    f= figure,
    hist(AAL(find(overlapF2>0)),unique(AAL(find(overlapF2>0))));
    % 30, 14, 16, (1), 34, 29, 68

    end
    
    %% training
    
    bodyk = (zscore(hearth.count_perday(keep,:))+zscore(adyad.count_perday(keep,:)))./2;
    bodyk = zscore(hearth.count_perday(keep,:))
    bt    = term(bodyk);
    
    M       = 1 + CH + A + S + bt + random(Sub) + I;
    slm     = SurfStatLinModS(Ck,M,SW);
    slm  	= SurfStatT(slm,(bodyk));
      
    f=figure,
    BoSurfStatViewData(slm.t,SInf,'')
    SurfStatColLim([-3 3])
    colormap([0 0 0; mycol.blackblue])
    %exportfigbo(f,[RPATH 'ATT_mod.png'],'png',10)
    clusatt    = BoSurfStatPlotAllStats(slm,SM, mask, 0.025, 0.05, [RPATH '/FIG1C.affect.hearth.cpd.WB' ],1);


  
end  

%%  ToM x PERSPECTIVE
for behavioralperspective = 1
    
    for figureC = 1
        keep1   = find(strcmp(groupN,'Perspective'))
        keep2   = intersect(find(X(:,1) >-666), find(sum(X,2)~=0));
        keep3   = find(abs(mean(X,2))<0.1);
        keep4   = find(~(strcmp(group4,'Group3')));
        mask    = mean(CT,1)>0.4;
        keep6   = find(tom>-666);
        keep6b  = find(abs(Tom)<666);
        keep    = mintersect(keep1,keep2,keep4,keep6,keep6b,keep3);
        size(mintersect(keep1,keep2,keep3, keep4)) % N = 120 (121)
        
        groupN(find(strcmp(groupN,'Control1'))) = {'Control'}
        
        size(keep) % N = 115
        size(mintersect(keep, find(strcmp(groupN,'Affect')))) %
        size(mintersect(keep, find(strcmp(groupN,'Perspective')))) % N = 115
        size(mintersect(keep, find(strcmp(groupN,'Presence'))))
        
        Ck      = X(keep,:);
        ik      = id(keep,:);
        ink     = idnum(keep);
        tk      = tpnum(keep,:);
        tnk     = tp(keep,:);
        gk      = group(keep);
        g4k     = group4(keep);
        gNk     = groupN(keep);
        ak      = age(keep);
        sk      = sex(keep);
        Tomk    = (Tom(keep));
        Tombfk  = (chtbl(keep));
        dk      = Days_last(keep);
        A       = term(ak);
        S       = term(sk);
        G       = term(gk);
        G4      = term(g4k);
        GN      = term(gNk);
        Sub     = term(ik);
        TO      = term(Tomk);
        T       = term(tk);
        Tn      = term(tnk);
        
        CHk      = mean(Ck,2);
        CH      = term(CHk);
        
        Hk    = Hear(keep);
        H     = term(Hk);
        
        Tbf  = term(Tombfk);
        
        D = term(dk);
        
        
        M       = 1 + CH + A + S + TO + random(Sub) + I;
        slm     = SurfStatLinModS(Ck,M,SW);
        slm  	= SurfStatT(slm,(Tomk));
        
        f=figure,
        BoSurfStatViewData(slm.t,SInf,'')
        SurfStatColLim([1.69 5])
        colormap([0 0 0; mycol.blackblue])
        exportfigbo(f,[RPATH 'FIG3C.TOM_mod.png'],'png',10)
        clustom   = BoSurfStatPlotAllStatsNew(slm,SM, mask, 0.025, 0.025, [RPATH '/FIG3C.all.tom.WB' ],1);
        clustom   = BoSurfStatPlotAllStats(slm,SInf, mask, 0.025, 0.025, [RPATH '/FIG3C.all.inf.tom.WB' ],1);
        
        %% test for complex model with random effects
        numclus = sum(clustom.clus.P < 0.025);
        for i = 1 : numclus
            seedx       = mean(CT(keep,clustom.clusid==i),2);
            seedxt      = term(seedx);
            
            M2       = 1 + CH + A + S + GN + TO + GN*TO + Sub;
            slm1 = SurfStatLinMod(seedx,M2);
            slm1 = SurfStatT(slm1, (GN.Perspective.*Tomk));
            disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
            % t = 2.56
            M3       = CH + A + S + (1 + TO)*(1 + GN)*(1 + random(Sub)) + I;
            slm1 = SurfStatLinMod(seedx,M3);
            slm1 = SurfStatT(slm1, (GN.Perspective.*Tomk));
            disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
            % t = 1.68
        end
        
        % AAL atlas
        numclus = sum(clustom.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(clustom.clusid==i);
            unique(AAL2C2)
            f= figure,
            hist(AAL(clustom.clusid==i),unique(AAL(clustom.clusid==i)));
            exportfigbo(f,[RPATH 'Tomxperspective' num2str(i) '.png'],'png',10);
            vertexid = find(clustom.t == max(clustom.t(clustom.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(clustom.t(clustom.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            
            %     1    57    59    61    67
            % t-val = 3.1114 AAL: 59, coordinate: -23.2231     -64.0538      48.9039
        end
        
        
        tomF3 = clustom.puncorr<0.025
        f = figure,
        BoSurfStatViewData(o_etombin.*tomF3, SInf,'')
        exportfigbo(f,[RPATH 'ROI.tom*tom.png'],'png',10);
        
        overlapF2 = o_etombin.*tomF3
        AAL2C2 = AAL(find(overlapF2>0));
        unique(AAL2C2)
        f= figure,
        hist(AAL(find(overlapF2>0)),unique(AAL(find(overlapF2>0))));
        % 64, 82, 68, 67
    end
    
    for  otheranalysis = 1
        %% training-hours
        
        bodyk = (zscore(thought.count_perday(keep,:)))+zscore(thought.count_perday(keep,:)))./2;
        bt    = term(bodyk);
        
        M       = 1 + CH + A + S + bt + random(Sub) + I;
        slm     = SurfStatLinModS(Ck,M,SW);
        slm  	= SurfStatT(slm,(bodyk));
        
        f=figure,
        BoSurfStatViewData(slm.t,SInf,'')
        SurfStatColLim([-3 3])
        colormap([0 0 0; mycol.blackblue])
        %exportfigbo(f,[RPATH 'ATT_mod.png'],'png',10)
        clusatt    = BoSurfStatPlotAllStats(slm,SM, mask, 0.025, 0.05, [RPATH '/FIG1C.perspective.thought.d7.WB' ],1);
        
        
        
        
        for j = [1 : numclus]
            target = mean(Ck(:,find(clustom.clusid==j)),2)
            M    = 1 + CH;
            slm = SurfStatLinMod(target,M);
            targetresid(keep,:) = target - slm.X*slm.coef;
            
            keepc1  = find(abs(Tom)<666);
            keepp = find(strcmp(groupN,'Perspective'));
            keepp = mintersect(keep,keepc1,keepp);
            f= figure,
            hold on
            scatter(Tom(keepp),targetresid(keepp),'Marker','o', 'MarkerFaceColor', [0 0.7 0], 'MarkerEdgeColor', [0 0 0]), lsline
            %axis([-0.3 0.3 -1 1])
            exportfigbo(f,[RPATH 'FIG3C.ToM.perspective.scatter' num2str(j) '.png'],'png',10)
            [r p] = corr(targetresid(keepp),Tom(keepp))
            % r = 0.31, p<0.0074
        end
        
        %% separate cohorts
        for j = [1 : numclus]
            target = mean(Ck(:,find(clustom.clusid==j)),2)
            M    = 1 + CH;
            slm = SurfStatLinMod(target,M);
            targetresid(keep,:) = target - slm.X*slm.coef;
            
            keepc1  = find(abs(Tom)<666);
            keepp1 = find(strcmp(groupN,'Perspective'));
            keepp2 = find(strcmp(group,'Group_1'));
            
            keepp = mintersect(keep,keepc1,keepp1,keepp2);
            f= figure,
            hold on
            scatter(Tom(keepp),targetresid(keepp),'Marker','o', 'MarkerFaceColor', [0 0.7 0], 'MarkerEdgeColor', [0 0 0]), lsline
            %axis([-0.3 0.3 -1 1])
            exportfigbo(f,[RPATH 'FIG3C.ToM.perspective.scatter.G1' num2str(j) '.png'],'png',10)
            [r p] = corr(targetresid(keepp),Tom(keepp))
            % r = 0.34, p<0.009
        end
        
        for j = [1 : numclus]
            target = mean(Ck(:,find(clustom.clusid==j)),2)
            M    = 1 + CH;
            slm = SurfStatLinMod(target,M);
            targetresid(keep,:) = target - slm.X*slm.coef;
            
            keepc1  = find(abs(Tom)<666);
            keepp1 = find(strcmp(groupN,'Perspective'));
            keepp2 = find(strcmp(group,'Group_2'));
            
            keepp = mintersect(keep,keepc1,keepp1,keepp2);
            f= figure,
            hold on
            scatter(Tom(keepp),targetresid(keepp),'Marker','o', 'MarkerFaceColor', [0 0.7 0], 'MarkerEdgeColor', [0 0 0]), lsline
            %axis([-0.3 0.3 -1 1])
            exportfigbo(f,[RPATH 'FIG3C.ToM.perspective.scatter.G2' num2str(j) '.png'],'png',10)
            [r p] = corr(targetresid(keepp),Tom(keepp))
            % r = 0.28, p<0.03
        end
        
        for i = 2:90
            for j = 2:90
                j
                myseed = mean(X(keep,AAL==i),2);
                mytarget = mean(X(keep,AAL==j),2);
                myseed(isnan(myseed))=1;
                mytarget(isnan(mytarget))=1;
                if i ~= j
                    Mb = 1  + CH + A + S + term(myseed);
                    slm = SurfStatLinMod(mytarget, Mb);
                    slm = SurfStatT(slm,myseed);
                    by_rsq = slm.df./(slm.t.^2)+1;
                    r = sqrt(1./by_rsq);
                    tval = slm.t;
                    test3(i-1,j-1) = tval;
                    test2(i-1,j-1) = r;
                else
                    test3(i-1,j-1) = 1;
                    test2(i-1,j-1) = 1;
                end
            end
        end
    end
end

%%  Make figure with the separate cohorts (Supplement?)
for FIG3C = 1
    keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
    keep1   = union(keep1,    find(strcmp(groupN,'Presence')));
    keep2   = intersect(find(X(:,1) >-666), find(sum(X,2)~=0));
    keep4   = find(~(strcmp(group4,'Group3')));
    mask    = mean(CT,1)>0.4;
    keep5   = find(tpnum>2);
    keep    = mintersect(keep1,keep2,keep4,keepMAIN,keep5);
    size(mintersect(keep1,keep2,keep3, keep4)) 

     
    size(keep) % 
    size(mintersect(keep, find(strcmp(groupN,'Affect')))) 
    size(mintersect(keep, find(strcmp(groupN,'Perspective')))) 
    size(mintersect(keep, find(strcmp(groupN,'Presence'))))

    Ck      = X(keep,:);
    ik      = id(keep,:);
    ink     = idnum(keep);
    tk      = tpnum(keep,:);
    tnk     = tp(keep,:);
    gk      = group(keep);
    g4k     = group4(keep);
    gNk     = groupN(keep);
    ak      = age(keep);
    sk      = sex(keep);

    dk      = Days_last(keep);
    A       = term(ak);
    S       = term(sk);
    G       = term(gk);
    G4      = term(g4k);
    GN      = term(gNk);
    Sub     = term(ik);

    T       = term(tk);
    Tn      = term(tnk);
    
    CHk      = mean(Ck,2);
    CH      = term(CHk);
   
    D = term(dk);
    
         
    %(in matlab_sofie/Project3) 
    for i = 1: 4
        target = mean(Ck(:,find( perspectiveF2.stats.clusid==i)),2);
        M      = 1 + CH;
        slm = SurfStatLinMod(target,M);
        test = target - slm.X*slm.coef;
        
         plotgroupsNEW(test, gNk, [RPATH 'SFx_scatterPerspectiveT2T3' num2str(i) ])
         
         M       = 1 + A + S + GN;
         slm     = SurfStatLinModS(test,M);
         slm  	 = SurfStatT(slm,(GN.Perspective-GN.Affect));
         slm.t
         
    end
    
     for i = 1: 4
        target = mean(Ck(:,find(affectiveF2.stats.clusid==i)),2);
        M      = 1 + CH;
        slm = SurfStatLinMod(target,M);
        test = target - slm.X*slm.coef;
        
         plotgroupsNEW(test, gNk, [RPATH 'SFx_scatterAffectT2T3' num2str(i) ])
         
         M       = 1 + A + S + GN;
         slm     = SurfStatLinModS(test,M);
         slm  	 = SurfStatT(slm,(GN.Affect-GN.Perspective));
         slm.t
         
    end
    
end

%% analysis that were done but not in current (Mai, 2016 version of paper). 

for FIG3A = 1

    keep1   = find(strcmp(groupN,'Presence'))
    keep4   = find(~(strcmp(group4,'Group3')));
    keep2   = intersect(find(X(:,1) >-100), find(sum(X,2)~=0));
    size(mintersect(keep1,keep2,keep4, find(strcmp(groupN,'Affect')))) % 
    size(mintersect(keep1,keep2,keep4, find(strcmp(groupN,'Perspective')))) 
    size(mintersect(keep1,keep2,keep4, find(strcmp(groupN,'Presence')))) % N = 133
    keep3   = find(abs(mean(X,2))<0.1);
    size(mintersect(keep1,keep2,keep3, keep4)) % N = 132

    keep5   = find(abs(Att)<666);
    keep6   = find(attention>-666)
    keep7   = find(abs(Att_c)<666);
    keep8   = find(attention_c>-666);
    keep    = mintersect(keep1,keep2,keep3,keep4,keep5);
 
    size(keep) % N = 102
    size(mintersect(keep, find(strcmp(groupN,'Affect')))) 
    size(mintersect(keep, find(strcmp(groupN,'Perspective')))) 
    size(mintersect(keep, find(strcmp(groupN,'Presence')))) % N = 102
    mean(Att(keep)) % 0.05
    std(Att(keep))  % 0.15
    mask    = mean(CT,1)>0.4;
    Ck      = X(keep,:);
    ik      = id(keep,:);
    ink     = idnum(keep);
    tk      = tpnum(keep,:);
    tnk     = tp(keep,:);
    dk      = Days_last(keep);
    gk      = group(keep);
    ak      = age(keep);
    sk      = sex(keep);
    Attk    = (Att(keep));
    Attblk  = (chAbl(keep));
    Att_ck  = (Att_c(keep));
    Attbl_ck= (attention_c(keep));
    A       = term(ak);
    S       = term(sk);
    G       = term(gk);
    G4      = term(g4k);
    %GN      = term(gNk);
    D       = term(dk);
    Sub     = term(ik);
    T       = term(tk);
    Tn      = term(tnk);
    ATT     = term(Attk);
    ATTBL   = term(Attblk);
    ATTc    = term(Att_ck);
    ATTBLc  = term(Attbl_ck);
    CTM     = mean(Ck,2);
    CH      = term(CTM);
    

    %% Figure 4A
    M       = 1 + CH + A + S + ATT  + random(Sub) + I;
    slm     = SurfStatLinModS(Ck,M,SW);
    slm  	= SurfStatT(slm,(Attk));
    
    f=figure,
    BoSurfStatViewData(slm.t,SInf,'')
    SurfStatColLim([-3 3])
    colormap([0 0 0; mycol.blackblue])
    exportfigbo(f,[RPATH 'ATT_mod.png'],'png',10)
    clusatt    = BoSurfStatPlotAllStatsNew(slm,SM, mask, 0.025, 0.025, [RPATH '/FIG1C.presence.att_CHk.WB' ],1);
    clusatt    = BoSurfStatPlotAllStats(slm,SInf, mask, 0.025, 0.025, [RPATH '/FIG1C.presence.inf.att_CHk.WB' ],1);
 
  
    %% body training
    bodyk = (zscore(body.count_perday(keep,:))) + zscore(body.count_perday(keep,:)))./2;
    bt    = term(bodyk);
    
    M       = 1 + CH + A + S + bt + random(Sub) + I;
    slm     = SurfStatLinModS(Ck,M,SW);
    slm  	= SurfStatT(slm,(bodyk));
      
    f=figure,
    BoSurfStatViewData(slm.t,SInf,'')
    SurfStatColLim([-3 3])
    colormap([0 0 0; mycol.blackblue])
    %exportfigbo(f,[RPATH 'ATT_mod.png'],'png',10)
    clusatt    = BoSurfStatPlotAllStats(slm,SM, mask, 0.025, 0.025, [RPATH '/FIG1C.presence.body.cpd.WB' ],1);


     
    numclus = sum(clusatt.clus.P < 0.025);
    for i = 1 : numclus
        seedx       = mean(X(keep,clusatt.clusid==i),2);
        seedxt      = term(seedx);
        
        M2       =1 + CH + A + S + GN + ATT + GN*ATT + Sub;
        slm1 = SurfStatLinMod(seedx,M2);
        slm1 = SurfStatT(slm1, (GN.Presence.*Attk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
% ROI: 0+1i,t= 2.6306
        M3       = CH + A + S + (1 + ATT)*(1 + GN)*(1 + random(Sub)) + I;
        slm1 = SurfStatLinModS(seedx,M3);
        slm1 = SurfStatT(slm1, (GN.Presence.*Attk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
% ROI: 0+1i,t= 2.0724
    end
    

    numclus = sum(clusatt.clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(clusatt.clusid==i);
        unique(AAL2C2)
        f= figure,
        hist(AAL(clusatt.clusid==i),unique(AAL(clusatt.clusid==i)));
        exportfigbo(f,[RPATH 'att-presence' num2str(i) '.png'],'png',10);
        vertexid = find(clusatt.t == max(clusatt.t(clusatt.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(clusatt.t(clusatt.clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
% NEW
%      1    55    81    85    89
% t-val = 3.5417 AAL: 89, coordinate: -58.2058     -32.3544     -18.8977
    end
    
    
    %% attention
    numclus = sum(clusatt.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(clusatt.clusid==j)),2);
         M    = 1 + CH; 
         slm = SurfStatLinMod(target,M); 
         targetresid = target - slm.X*slm.coef; 
         f= figure,
         hold on
         scatter(Attk,target, 60, 'Marker','o', 'MarkerFaceColor', [0.9 0.9 0], 'MarkerEdgeColor', [0 0 0]), lsline
         %axis([-0.3 0.3 -1 1])
         exportfigbo(f,[RPATH 'FIG3C.Attadj.Presence.scatter' num2str(j) '.png'],'png',10)
         [r p] = corr(Attk,target.*Attblk)
%          r = 0.34
    end
    
    %% separate groups 
    numclus = sum(clusatt.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(clusatt.clusid==j)),2);
         M    = 1 + CH + A + S; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepatt1  = find(abs(Att)<666);
         keepp = find(strcmp(groupN,'Presence'));
         keepp2 = find(strcmp(group,'Group_1'));
         keepatt = mintersect(keep,keepatt1,keepp,keepp2);
         f= figure,
         hold on
         scatter(Att(keepatt),targetresid(keepatt), 60, 'Marker','o', 'MarkerFaceColor', [0.9 0.9 0], 'MarkerEdgeColor', [0 0 0]), lsline
         axis([-0.3 0.6 -0.2 0.2])
         exportfigbo(f,[RPATH 'FIG3C.Attch.Presence.scatter.G1.' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepatt),Att(keepatt))
%          r = 0.45   p =   0.0007
    end
    
    numclus = sum(clusatt.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(clusatt.clusid==j)),2);
         M    = 1 + CH + A + S; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepatt1  = find(abs(Att)<666);
         keepp = find(strcmp(groupN,'Presence'));
         keepp2 = find(strcmp(group,'Group_2'));
         keepatt = mintersect(keep,keepatt1,keepp,keepp2);
         f= figure,
         hold on
         scatter(Att(keepatt),targetresid(keepatt), 60, 'Marker','o', 'MarkerFaceColor', [0.9 0.9 0], 'MarkerEdgeColor', [0 0 0]), lsline
         axis([-0.3 0.6 -0.2 0.2])
         exportfigbo(f,[RPATH 'FIG3C.Attch.Presence.scatter.G2.' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepatt),Att(keepatt))
%          r = 0.17
    end
    

end

%% quick graph with aparc - aal 
for FIG3C = 1
    keep1   = union(find(strcmp(groupN,'Affect')),    find(strcmp(groupN,'Perspective')));
    keep1   = union(keep1,    find(strcmp(groupN,'Presence')));
    keep2   = intersect(find(X(:,1) >-666), find(sum(X,2)~=0));
    keep3   = find(abs(mean(X,2))<0.1);
    keep4   = find(~(strcmp(group4,'Group3')));
    mask    = mean(CT,1)>0.4;

    keep    = mintersect(keep1,keep2,keep4,keep3);
    size(mintersect(keep1,keep2,keep3, keep4)) % N = 120 (121)

    groupN(find(strcmp(groupN,'Control1'))) = {'Control'}
     
    size(keep) % N = 115
    size(mintersect(keep, find(strcmp(groupN,'Affect')))) % 
    size(mintersect(keep, find(strcmp(groupN,'Perspective')))) % N = 115
    size(mintersect(keep, find(strcmp(groupN,'Presence'))))

    Ck      = X(keep,:);
    ik      = id(keep,:);
    ink     = idnum(keep);
    tk      = tpnum(keep,:);
    tnk     = tp(keep,:);
    gk      = group(keep);
    g4k     = group4(keep);
    gNk     = groupN(keep);
    ak      = age(keep);
    sk      = sex(keep);

    dk      = Days_last(keep);
    A       = term(ak);
    S       = term(sk);
    G       = term(gk);
    G4      = term(g4k);
    GN      = term(gNk);
    Sub     = term(ik);

    T       = term(tk);
    Tn      = term(tnk);
    
    CHk      = mean(Ck,2);
    CH      = term(CHk);
    
   
    
    D = term(dk);
    

 
    for i = 1:90
        for j = 1:90
            j
            myseed = mean(X(keep,AAL==i),2);
            mytarget = mean(X(keep,AAL==j),2);
            myseed(isnan(myseed))=1;
            mytarget(isnan(mytarget))=1;
            if i ~= j
                Mb = 1  + CH + A + S + GN + term(myseed) + GN*term(myseed);
                slm = SurfStatLinMod(mytarget, Mb);
                slm = SurfStatT(slm,((GN.Affect.*myseed)-(GN.Perspective.*myseed))-(GN.Presence.*myseed));
                by_rsq = slm.df./(slm.t.^2)+1;
                r = sqrt(1./by_rsq);
                tval = slm.t;
                test3(i,j) = tval;
                test2(i,j) = r;
            else
                test3(i,j) = 1;
                test2(i,j) = 1;
            end
        end
    end
    %% binarize test3
    test4 = zeros(90,90);
    test4(test3>2) = 1;
    test4l = sum(test4);
    for i = 1:20484
        i
        if AAL(i) ~= 0
        vertexnew(i,:) = test4l(AAL(i));
        else
            vertexnew(i,:) = 0;
        end
    end
    f =figure,
    BoSurfStatViewData(vertexnew,SInf,'')


  for i = 2:90
        for j = 2:90
            j
            myseed = mean(X(keep,AAL==i),2);
            mytarget = mean(X(keep,AAL==j),2);
            myseed(isnan(myseed))=1;
            mytarget(isnan(mytarget))=1;
            if i ~= j 
            Mb = 1  + CH + A + S + term(myseed);
            slm = SurfStatLinMod(mytarget, Mb);
            slm = SurfStatT(slm,myseed.*Compk);
            by_rsq = slm.df./(slm.t.^2)+1;
            r = sqrt(1./by_rsq);
            tval = slm.t;
            test3(i-1,j-1) = tval;
            test2(i-1,j-1) = r;
            else
                test3(i-1,j-1) = 1;
                test2(i-1,j-1) = 1;
            end
        end
    end
end


%% TC3 change NOT USED

for FIG12 = 1
    keep    = mintersect(find(isthere_ct), ...
        find(tp_t0_there==1), ...
        find(tp_t1_there==1), ...
        union(find(tpnum==0), find(tpnum==1)));
    
    size(keep) % N = 716
    mask    = mean(CT,1)>0.4;
    Ck      = CT(keep,:);
    ik      = id(keep,:);
    ink     = idnum(keep);
    tk      = tpnum(keep,:);
    dk      = daysfrombl(keep);
    dkfac   = tp(keep);
    gk      = group(keep);
    g4k     = group4(keep);
    ak      = age(keep);
    sk      = sex(keep);
    
    A       = term(ak);
    S       = term(sk);
    G       = term(gk);
    G4      = term(g4k);
    D       = term(dk);
    Sub     = term(ik);
    DT      = term(dkfac);  
    
    CHbl    = CT(keep,:);
    CHk =  mean(CHbl,2);
    CH  =  term(CHk);
    % model2
    M       = 1  + CH + A + S + G4 + DT + G4*DT + random(Sub) + I;
    slm = SurfStatLinModS(Ck,M,SW);
    
    %% perspective - affect controlled for presence
    contrast = ((G4.Group3.* DT.T1)-(G4.Group3.* DT.T0)) - ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) ...
        - ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0))
    % FIGURE2A -- TC2 - TC2
    slm  	= SurfStatT(slm,contrast);
    perspectiveF2.stats   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.TC3-TC12.FAC.days.M3.'],1);
    perspectiveF2.stats   = BoSurfStatPlotAllStats(slm, SInf, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.inf.TC3-TC12.FAC.days.'],1);
    
    %% test for complex model with random effects // FIXED effects
    numclus = sum(perspectiveF2.stats.clus.P < 0.025);
    for i = 1 : numclus
        seedx       = mean(CT(keep,perspectiveF2.stats.clusid==i),2);
        seedxt      = term(seedx);
        
        M2       = 1 + CH + A + S + G4 + DT + G4*DT + Sub;
        slm1 = SurfStatLinModS(seedx,M2);
        slm1 = SurfStatT(slm1, contrast);
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
% ROI: 1,t= 4.0308
% ROI: 1,t= 4.4525
% ROI: 1,t= 3.7377
% ROI: 1,t= 3.5254
%         M3       = CH + A + S + (1 + G4)*(1 + DT)*(1 + random(Sub)) + I;
%         slm1 = SurfStatLinMod(seedx,M3);
%         slm1 = SurfStatT(slm1, contrast);
%         disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
    end
    
    %% AAL
%     AAL2C2 = AAL(perspectiveF2.stats.clusid==3);
%     unique(AAL2C2)
%     f= figure,
%     hist(AAL(perspectiveF2.stats.clusid==1),unique(AAL(perspectiveF2.stats.clusid==1)))
%     exportfigbo(f,[RPATH 'F2BperspectiveAAL' num2str(i) '.png'],'png',10)
    
    numclus = sum(perspectiveF2.stats.clus.P < 0.025);
    for i = 1 : numclus
        vertexid = find(perspectiveF2.stats.t == max(perspectiveF2.stats.t(perspectiveF2.stats.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(perspectiveF2.stats.t(perspectiveF2.stats.clusid==i));
        AAL2C2 = AAL(perspectiveF2.stats.clusid==i);
        unique(AAL2C2)
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         0    81    85    89
%         t-val = 3.6348 AAL: 85, coordinate: -52.8215     -48.3212      3.25627
%         0    43    44    46    48    50
%         t-val = 3.1379 AAL: 48, coordinate: 10.4429     -84.4072     -13.6755
%         0    43    45    47    49    51
%         t-val = 3.0483 AAL: 43, coordinate: -10.8362     -86.9299      7.04825

%NEW
%      1    81    85    89
% t-val = 3.6375 AAL: 85, coordinate: -51.2848     -43.3906     0.232559
%      1    82    86    90
% t-val = 4.0353 AAL: 86, coordinate: 65.302     -42.8761     -8.28508 
%      1    43    44    46    48    50
% t-val = 3.0493 AAL: 48, coordinate: 10.4429     -84.4072     -13.6755
%      1    43    45    47    49    51
% t-val = 3.0419 AAL: 43, coordinate: -10.8362     -86.9299      7.04825
    end
      
    % cluster significance
    % cohens D:
    numclus = sum(perspectiveF2.stats.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,perspectiveF2.stats.clusid==i),2);
        slm      = SurfStatLinMod(ROI,M);
        slm      = SurfStatT(slm,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slm.t)])   
        d = 2*(slm.t) / sqrt(slm.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(perspectiveF2.stats.clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val 3.9723
% Cohen d 0.30006
% cluster sig 3.0145e-05
% *******************************************************
% Post-hoc Clus: 2
% t-val 4.3946
% Cohen d 0.33196
% cluster sig 0.00032642
% *******************************************************
% Post-hoc Clus: 3
% t-val 3.7023
% Cohen d 0.27967
% cluster sig 0.00039396
% *******************************************************
% Post-hoc Clus: 4
% t-val 3.515
% Cohen d 0.26552
% cluster sig 0.0011375
    end
    
     %% affect -perspective controlled for presence
    contrast = (((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) + ((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2))) - ...
        (((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) + ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2))) - ...
            (((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_2.* DT.T0)))
    slm  	= SurfStatT(slm,contrast);
    affectiveF2.stats   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.affect-perspective.b.FAC.M3.days.'],1);
    affectiveF2.stats   = BoSurfStatPlotAllStats(slm, SInf, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.INF.affect-perspective.b.FAC.days.'],1);
   
    %% test for complex model with random effects
    numclus = sum(affectiveF2.stats.clus.P < 0.025);
    for i = 1 : numclus
        seedx       = mean(CT(keep,affectiveF2.stats.clusid==i),2);
        seedxt      = term(seedx);
        
        M2       = 1 + CH + A + S + G4 + DT + G4*DT + Sub;
        slm1 = SurfStatLinMod(seedx,M2);
        slm1 = SurfStatT(slm1, contrast);
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
% ROI: 1,t= 5.6359
% ROI: 1,t= 4.1286
% ROI: 1,t= 3.7002
% ROI: 1,t= 3.1319
    end
    
    %% AAL
    numclus = sum(affectiveF2.stats.clus.P < 0.025);
    for i = 1 : numclus
        vertexid = find(affectiveF2.stats.t == max(affectiveF2.stats.t(affectiveF2.stats.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(affectiveF2.stats.t(affectiveF2.stats.clusid==i));
        AAL2C2 = AAL(affectiveF2.stats.clusid==i);
        unique(AAL2C2)
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         2     8    12    18    58
%         t-val = 3.9387 AAL: 2, coordinate: 44.9294       1.0018      46.7104
%         8    12    14    16    18
%         t-val = 4.1148 AAL: 14, coordinate: 54.2935       29.129      2.25135
%         64    82
%         t-val = 2.8816 AAL: 82, coordinate: 64.0045     -33.1973      16.0255
%         0    33    35    67    69
%         t-val = 3.4546 AAL: 33, coordinate: -7.97467     -40.1174      38.1422
%         0    61    63    65    85
%         t-val = 3.2721 AAL: 61, coordinate: -51.5047     -50.0944      36.4772

%NEW
%      1     2     8    12    14    16    18    30    58    64
% t-val = 4.1065 AAL: 14, coordinate: 53.3627      27.3439       1.9629
%      1    33    35    67    69
% t-val = 3.4612 AAL: 33, coordinate: -7.97467     -40.1174      38.1422
%      1    61    63    65    85
% t-val = 3.3241 AAL: 61, coordinate: -51.5047     -50.0944      36.4772
%      1    51    53    55    85    89
% t-val = 2.9988 AAL: 85, coordinate: -63.1341     -38.1133      -14.029
    end
    

      
    % cluster significance
    % cohens D:
    numclus = sum(affectiveF2.stats.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,affectiveF2.stats.clusid==i),2);
        slm1      = SurfStatLinMod(ROI,M);
        slm1      = SurfStatT(slm1,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slm1.t)])   
        d = 2*(slm1.t) / sqrt(slm1.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(affectiveF2.stats.clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val 5.5822
% Cohen d 0.42167
% cluster sig 1.2376e-06
% *******************************************************
% Post-hoc Clus: 2
% t-val 4.1021
% Cohen d 0.30987
% cluster sig 0.0033708
% *******************************************************
% Post-hoc Clus: 3
% t-val 3.632
% Cohen d 0.27436
% cluster sig 0.013331
% *******************************************************
% Post-hoc Clus: 4
% t-val 3.1256
% Cohen d 0.2361
% cluster sig 0.018762
    end
    
    
    %% To-T2 
    contrast = ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T0)) - ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T0)) 
    slm  	= SurfStatT(slm,contrast);
    T2T0stats   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.T2-T0.affect-perspective.FAC.days.'],1);
     T2T0stats   = BoSurfStatPlotAllStatsNew(slm, SInf, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.INF.T2-T0.affect-perspective.FAC.days.'],1);
    F2C_stats = stats;
    
     %% test for complex model with random effects
    numclus = sum(T2T0stats.clus.P < 0.025);
    for i = 1 : numclus
        seedx       = mean(CT(keep,T2T0stats.clusid==i),2);
        seedxt      = term(seedx);
        
        M2       = 1 + CH + A + S + G4 + DT + G4*DT + Sub;
        slm1 = SurfStatLinMod(seedx,M2);
        slm1 = SurfStatT(slm1, contrast);
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
% ROI: 1,t= 3.9804
% ROI: 1,t= 4.2322
    end
    
    %% test for complex model with random effects
    numclus = sum(T2T0stats.inv_clus.P < 0.025);
    for i = 1 : numclus
        seedx       = mean(CT(keep,T2T0stats.inv_clusid==i),2);
        seedxt      = term(seedx);
        
        M2       = 1 + CH + A + S + G4 + DT + G4*DT + Sub;
        slm1 = SurfStatLinMod(seedx,M2);
        slm1 = SurfStatT(slm1, contrast);
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
% ROI: 1,t= -4.7594
% ROI: 1,t= -4.7681
    end
    
    
    %% Affect - perspective
    numclus = sum(T2T0stats.clus.P < 0.025);
    for i = 1 : numclus
        vertexid = find(T2T0stats.t == max(T2T0stats.t(T2T0stats.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(T2T0stats.t(T2T0stats.clusid==i));
        AAL2C2 = AAL(T2T0stats.clusid==i);
        unique(AAL2C2)
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         0     3     7    13    19    23    31    33
%         t-val = 3.7661 AAL: 23, coordinate: -12.0624      36.2348      22.4073
%         0    19    33    67    69
%         t-val = 3.4457 AAL: 33, coordinate: -11.2437     -3.47004       43.516
%         0    12    14    18    30
%         t-val = 3.6728 AAL: 12, coordinate: 41.2897      5.83642      8.57816

%NEW

%      1     3     7    13    19    23    31    33
% t-val = 3.8197 AAL: 23, coordinate: -12.0624      36.2348      22.4073
%      1    19    33    67    69
% t-val = 3.4648 AAL: 33, coordinate: -11.2437     -3.47004       43.516
    end
      
    % cluster significance
    % cohens D:
    numclus = sum(T2T0stats.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,T2T0stats.clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(T2T0stats.clus.P(i))])
%         *******************************************************
%         Post-hoc Clus: 1
%         t-val 4.012
%         Cohen d 0.30306
%         cluster sig 4.3115e-06
%         *******************************************************
%         Post-hoc Clus: 2
%         t-val 4.2546
%         Cohen d 0.32139
%         cluster sig 0.006329
    end
    
    %% Perspective - affect
    numclus = sum(T2T0stats.inv_clus.P < 0.025);
    for i = 1 : numclus
        vertexid = find(T2T0stats.inv_t == max(T2T0stats.inv_t(T2T0stats.inv_clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(T2T0stats.inv_t(T2T0stats.inv_clusid==i));
        AAL2C2 = AAL(T2T0stats.inv_clusid==i);
        unique(AAL2C2)
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         0    17    57    61    63    65    81    85    89
%         t-val = 4.1983 AAL: 85, coordinate: -52.8215     -48.3212      3.25627
%         0    82    84    86    88    90
%         t-val = 3.9859 AAL: 82, coordinate: 52.2421     -7.97336     -15.1101
% NEW

% 1    17    57    61    63    65    81    85    89
% t-val = 4.2654 AAL: 81, coordinate: -64.3335     -39.6449        13.92
% 1    82    84    86    88    90
% t-val = 4.9318 AAL: 86, coordinate: 56.1768     -16.6383     -15.5999
    end
      
    % cluster significance
    % cohens D:
    numclus = sum(T2T0stats.inv_clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,T2T0stats.inv_clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(T2T0stats.inv_clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val -4.7749
% Cohen d -0.36069
% cluster sig 1.2528e-06
% *******************************************************
% Post-hoc Clus: 2
% t-val -4.7491
% Cohen d -0.35874
% cluster sig 4.8931e-05
    end
    
    
    %% PRESENCE G2
    
    contrast = ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) - ((G4.Control.* DT.T1)-(G4.Control.* DT.T0))))    
    slm  	= SurfStatT(slm,contrast);
    statsG2   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.G2.presence-control.T0T1.FAC.days.'],1);
   
        %% post-hoc 
    % cluster significance
    % cohens D:
    numclus = sum(statsG2.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statsG2.clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statsG2.clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val 5.258
% Cohen d 0.39718
% cluster sig 1.2528e-06
    end
    
    numclus = sum(statsG2.inv_clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statsG2.inv_clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statsG2.inv_clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val -4.981
% Cohen d -0.37626
% cluster sig 1.2667e-06
% *******************************************************
% Post-hoc Clus: 2
% t-val -3.6343
% Cohen d -0.27453
% cluster sig 0.0050648
% *******************************************************
% Post-hoc Clus: 3
% t-val -3.5126
% Cohen d -0.26534
% cluster sig 0.009902
    end
    
    numclus = sum(statsG2.clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statsG2.clusid==i);
        unique(AAL2C2)
        %f= figure,
        %hist(AAL(statsG2.clusid==i),unique(AAL(statsG2.clusid==i)));
        %exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
        vertexid = find(statsG2.t == max(statsG2.t(statsG2.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statsG2.t(statsG2.clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)]) 
%         0     4     6     8    10    14    16    20    24    26    28    32    34
%         t-val = 4.659 AAL: 16, coordinate: 44.6538      48.6399     -7.64978

%NEW
%      1     4     6     8    10    14    16    20    24    26    28    32    34
% t-val = 4.6543 AAL: 16, coordinate: 44.6538      48.6399     -7.64978
    end
    
    numclus = sum(statsG2.inv_clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statsG2.inv_clusid==i);
        unique(AAL2C2)
        %f= figure,
        %hist(AAL(statsG2.inv_clusid==i),unique(AAL(statsG2.inv_clusid==i)));
        %exportfigbo(f,[RPATH 'G2presenceN' num2str(i) '.png'],'png',10);
        vertexid = find(statsG2.inv_t == max(statsG2.inv_t(statsG2.inv_clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statsG2.inv_t(statsG2.inv_clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         40    44    46    48    68
%         t-val = 4.3499 AAL: 44, coordinate: 15.6613     -56.3666      8.03368
%         0    18    30    64    82
%         t-val = 3.4321 AAL: 64, coordinate: 59.5053     -36.5102      23.0216
%         7     9    13    15
%         t-val = 3.545 AAL: 13, coordinate: -40.1846      29.9264      13.0839

%NEW
%     40    44    46    48    68
% t-val = 4.3292 AAL: 44, coordinate: 15.6613     -56.3666      8.03368
%      1    18    30    64    80    82
% t-val = 3.4476 AAL: 64, coordinate: 59.5053     -36.5102      23.0216
%      7     9    13    15
% t-val = 3.571 AAL: 13, coordinate: -40.1846      29.9264      13.0839

    end
    
    contrast = ((((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) - ((G4.Control.* DT.T1)-(G4.Control.* DT.T0)))) 
    
    slm  	= SurfStatT(slm,contrast);
    statsG1   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.G1.presence-control.T0T1.FAC.days.'],1);
    
    %% post-hoc 
    % cluster significance
    % cohens D:
    numclus = sum(statsG1.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statsG1.clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statsG1.clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val 6.2705
% Cohen d 0.47367
% cluster sig 1.2528e-06
% *******************************************************
% Post-hoc Clus: 2
% t-val 6.3206
% Cohen d 0.47745
% cluster sig 1.2528e-06
    end
    
    numclus = sum(statsG1.inv_clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statsG1.inv_clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statsG1.inv_clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val -6.134
% Cohen d -0.46335
% cluster sig 1.2528e-06
% *******************************************************
% Post-hoc Clus: 2
% t-val -5.363
% Cohen d -0.40511
% cluster sig 1.2528e-06
% *******************************************************
% Post-hoc Clus: 3
% t-val -5.6242
% Cohen d -0.42484
% cluster sig 1.2529e-06
% *******************************************************
% Post-hoc Clus: 4
% t-val -4.406
% Cohen d -0.33283
% cluster sig 5.5894e-05
% *******************************************************
% Post-hoc Clus: 5
% t-val -4.066
% Cohen d -0.30714
% cluster sig 0.020707
    end
    
    numclus = sum(statsG1.clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statsG1.clusid==i);
        unique(AAL2C2)
        vertexid = find(statsG1.t == max(statsG1.t(statsG1.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statsG1.t(statsG1.clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         0     3     5     7     9    11    13    15    19    23    25    27    29    31    33
%         t-val = 4.9987 AAL: 3, coordinate: -15.1788      49.6029       36.329
%         0     4     6     8    10    16    20    24    26    28    32    34
%         t-val = 4.9618 AAL: 24, coordinate: 6.19002      27.0555      58.8839

%NEW

%      1     3     5     7     9    11    13    15    19    23    25    27    29    31    33
% t-val = 5.0196 AAL: 3, coordinate: -15.1788      49.6029       36.329
%      1     4     6     8    10    16    20    24    26    28    32    34    68
% t-val = 5.072 AAL: 24, coordinate: 6.19002      27.0555      58.8839
    end
    
    numclus = sum(statsG1.inv_clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statsG1.inv_clusid==i);
        unique(AAL2C2)
        f= figure,
        hist(AAL(statsG1.inv_clusid==i),unique(AAL(statsG1.inv_clusid==i)));
        exportfigbo(f,[RPATH 'G2presenceN' num2str(i) '.png'],'png',10);
        vertexid = find(statsG1.inv_t == max(statsG1.inv_t(statsG1.inv_clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statsG1.inv_t(statsG1.inv_clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         0    51    53    55    57    61    63    65    81    85    89
%         t-val = 5.3851 AAL: 85, coordinate: -56.267     -45.6577     -12.1994
%         0    82    86    88    90
%         t-val = 4.0975 AAL: 86, coordinate: 60.0739     0.359326     -26.7405
%         0    40    43    44    46    48    50    52    54    68
%         t-val = 4.0536 AAL: 44, coordinate: 7.81859     -96.7332      7.14824
%         2    12    14    18    58
%         t-val = 3.3807 AAL: 18, coordinate: 60.3209      3.89001      12.8612
%         64    66    82
%         t-val = 3.1803 AAL: 64, coordinate: 62.7818     -41.4796       22.412

%NEW
%      1    51    53    55    57    61    63    65    81    85    89
% t-val = 5.3919 AAL: 85, coordinate: -56.267     -45.6577     -12.1994
%      1    40    43    44    46    48    50    52    54    68
% t-val = 4.057 AAL: 44, coordinate: 7.81859     -96.7332      7.14824
%      1    64    66    82    84    86    88    90
% t-val = 4.3105 AAL: 86, coordinate: 57.2689     -2.11664     -31.0907
%      2    12    14    18    58    64
% t-val = 3.6607 AAL: 64, coordinate: 58.8944     -21.1868      27.0442
%      1    43    45    47    49    51    53
% t-val = 4.1444 AAL: 51, coordinate: -16.61529     -100.7017       4.65367
    end
    
    
    contrast = ((((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) + ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0))) ...
        - ((G4.Control.* DT.T1)-(G4.Control.* DT.T0))) 
    
    slm  	= SurfStatT(slm,contrast);
    statsG1   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.presence-control.T0T1.FAC.days.'],1);
  
    
    
    
    inter_tval = (statsG1.t) .* (statsG2.t)
    inter_len = (statsG1.p<p_thresh) .* (statsG2.p<p_thresh)
    f=figure;
    BoSurfStatViewData(double(inter_len>0),SInf,'');
    colormap([0 0 0; 1 0 0])
    exportfigbo(f,[RPATH 'FIGURE1.intersect_TC1-RCC_and_TC2-RCC_len.png'],'png',10);
    
    inter_neg_len = (statsG1.inv_p<p_thresh) .* (statsG2.inv_p<p_thresh)
    f=figure;
    BoSurfStatViewData(double(inter_neg_len>0),SInf,'');
    colormap([0 0 0; 0 0 1])
    exportfigbo(f,[RPATH 'FIGURE1.intersect_TC1-RCC_and_TC2-RCC_invlen.png'],'png',10);
   
    % FIGURE1F -- OVERLAP
    inter_str = (statsG1.pval.C<pwe) .* (statsG2.pval.C<pwe)
    f=figure;
    BoSurfStatViewData(double(inter_str>0),SM,'');
    colormap([0 0 0; 1 0 0])
    exportfigbo(f,[RPATH 'FIGURE1.intersect_TC1-RCC_and_TC2-RCC_str.png'],'png',10);
    
    inter_neg_str = (statsG1.inv_pval.C<pwe) .* (statsG2.inv_pval.C<pwe)
    f=figure;
    BoSurfStatViewData(double(inter_neg_str>0),SM,'');
    colormap([0 0 0; 0 0 1])
    exportfigbo(f,[RPATH 'FIGURE1.intersect_TC1-RCC_and_TC2-RCC_invstr.png'],'png',10);
   
    %% S1: T1_T2: 
    
    contrast = ((((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) - ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)))) - ...
        ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) - ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))
    
    slm  	= SurfStatT(slm,contrast);
    statsG2G1_T21   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.G2-G1.T1T2+.FAC.days.'],1);
  
    numclus = sum(statsG2G1_T21.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statsG2G1_T21.clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statsG2G1_T21.clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val 5.2439
% Cohen d 0.39612
% cluster sig 1.2547e-06
    end
    
    numclus = sum(statsG2G1_T21.clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statsG2G1_T21.clusid==i);
        unique(AAL2C2)
        %f= figure,
        %hist(AAL(statsG2G1_T21.clusid==i),unique(AAL(statsG2G1_T21.clusid==i)));
        %exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
        vertexid = find(statsG2G1_T21.t == max(statsG2G1_T21.t(statsG2G1_T21.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statsG2G1_T21.t(statsG2G1_T21.clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])    
        %  0     1     3     5     7     9    11    13    15    23    25
        % t-val = 5,0324 AAL: 13, coordinate: -40.1846      29.9264      13.0839
        % NEW
%         1     3     5     7     9    11    13    15    23    25
%         t-val = 5.0486 AAL: 13, coordinate: -40.1846      29.9264      13.0839
    end
    
    
    contrast = ((((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1))))- ...
        (((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)))
    
    slm  	= SurfStatT(slm,contrast);
    statsG1G2_T21   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.G1-G2.T1T2+.FAC.days.'],1);
    
    numclus = sum(statsG1G2_T21.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statsG1G2_T21.clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statsG1G2_T21.clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val 5.7655
% Cohen d 0.43552
% cluster sig 1.2527e-06
% *******************************************************
% Post-hoc Clus: 2
% t-val 3.6182
% Cohen d 0.27331
% cluster sig 8.8456e-05
% *******************************************************
% Post-hoc Clus: 3
% t-val 3.9276
% Cohen d 0.29669
% cluster sig 0.01308
% *******************************************************
% Post-hoc Clus: 4
% t-val 3.6801
% Cohen d 0.27799
% cluster sig 0.017573
    end
    
    numclus = sum(statsG1G2_T21.clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statsG1G2_T21.clusid==i);
        unique(AAL2C2)
%         f= figure,
%         hist(AAL(statsG1G2_T21.clusid==i),unique(AAL(statsG1G2_T21.clusid==i)));
%         exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
        vertexid = find(statsG1G2_T21.t == max(statsG1G2_T21.t(statsG1G2_T21.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statsG1G2_T21.t(statsG1G2_T21.clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])    
%         2     8    12    14    16    18    58
%         t-val = 4.10 AAL: 14, coordinate: 54.2935       29.129      2.25135
%         64    82
%         t-val = 3.17 AAL: 82, coordinate: 64.0045     -33.1973      16.0255
%         0    51    53    55    85    89
%         t-val = 3.61 AAL: 85, coordinate: -58.8292     -44.6253     -10.4067
%         0    61    63    65    85
%         t-val = 3.20 AAL: 61, coordinate: -57.5676     -44.5012      39.0814
%         0    33    35    67    69
%         t-val = 3.57 AAL: 33, coordinate: -7.97467     -40.1174      38.1422

%NEW
%      1     2     8    12    14    16    18    30    58    64
% t-val = 4.09 AAL: 14, coordinate: 54.2935       29.129      2.25135
%      1    51    53    55    85    89
% t-val = 3.6237 AAL: 85, coordinate: -58.8292     -44.6253     -10.4067
%      1    33    35    67    69
% t-val = 3.5742 AAL: 33, coordinate: -7.97467     -40.1174      38.1422
%      1    61    63    65    85
% t-val = 3.264 AAL: 61, coordinate: -57.5676     -44.5012      39.0814
    end
    
    
    %% T3_T2: 
    
    contrast = ((((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)) - ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)))) - ...
        ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1))- ...
        ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) - ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))
    
    slm  	= SurfStatT(slm,contrast);
    statsG2G1_T23   = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.G2-G1.T3T2+.FAC.days.'],1);
    
    
    numclus = sum(statsG2G1_T23.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statsG2G1_T23.clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statsG2G1_T23.clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val 3.6585
% Cohen d 0.27636
% cluster sig 0.00085688
% *******************************************************
% Post-hoc Clus: 2
% t-val 4.6712
% Cohen d 0.35285
% cluster sig 0.0011961
% *******************************************************
% Post-hoc Clus: 3
% t-val 3.9191
% Cohen d 0.29604
% cluster sig 0.010455
% *******************************************************
% Post-hoc Clus: 4
% t-val 4.2289
% Cohen d 0.31945
% cluster sig 0.013548
    end
    
    numclus = sum(statsG2G1_T23.clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statsG2G1_T23.clusid==i);
        unique(AAL2C2)
%         f= figure,
%         hist(AAL(statsG2G1_T23.clusid==i),unique(AAL(statsG2G1_T23.clusid==i)));
%         exportfigbo(f,[RPATH 'G2_G1_T2T3' num2str(i) '.png'],'png',10);
        vertexid = find(statsG2G1_T23.t == max(statsG2G1_T23.t(statsG2G1_T23.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statsG2G1_T23.t(statsG2G1_T23.clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])    
%      0    18    62    64    66    82    86
% t-val = 3.5695 AAL: 82, coordinate: 59.4246     -43.6284      16.3785
%     0    19    33    67    69
% t-val = 4.3564 AAL: 33, coordinate: -17.5391     -33.6701       41.906
%     40    56
% t-val = 3.4905 AAL: 40, coordinate: 23.7569     -29.2995     -22.9605
%     39    47    55
% t-val = 3.9825 AAL: 39, coordinate: -22.6042     -31.1043     -22.4635

%NEW

%   64    66    82    86
% t-val = 3.3937 AAL: 66, coordinate: 47.9344     -54.9842      31.4392
%      1    19    33    67    69
% t-val = 4.417 AAL: 33, coordinate: -17.5391     -33.6701       41.906
%     39    47    55
% t-val = 4.0056 AAL: 39, coordinate: -22.6042     -31.1043     -22.4635
%     40    56
% t-val = 3.3402 AAL: 40, coordinate: 21.515     -29.2189     -20.7169
    end
    
    
    contrast = ((((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)) - ((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)))) - ...
        ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) - ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1))- ...
        ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) - ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))
    
    slm  	= SurfStatT(slm,contrast);
    statsG1G2_T23  = BoSurfStatPlotAllStats(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.G1-G2.T3T2+.FAC.days.'],1);
    
    numclus = sum(statsG1G2_T23.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statsG1G2_T23.clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statsG1G2_T23.clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val 4.4961
% Cohen d 0.33963
% cluster sig 1.305e-06
% *******************************************************
% Post-hoc Clus: 2
% t-val 3.774
% Cohen d 0.28508
% cluster sig 7.1078e-05
% *******************************************************
% Post-hoc Clus: 3
% t-val 4.0503
% Cohen d 0.30595
% cluster sig 0.00010615
% *******************************************************
% Post-hoc Clus: 4
% t-val 3.846
% Cohen d 0.29052
% cluster sig 0.0026042
    end
    
    numclus = sum(statsG1G2_T23.clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statsG1G2_T23.clusid==i);
        unique(AAL2C2)
        f= figure,
        hist(AAL(statsG1G2_T23.clusid==i),unique(AAL(statsG1G2_T23.clusid==i)));
        exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
        vertexid = find(statsG1G2_T23.t == max(statsG1G2_T23.t(statsG1G2_T23.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statsG1G2_T23.t(statsG1G2_T23.clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)]) 
%         0    43    45    47    49    55
%         t-val = 3.8638 AAL: 43, coordinate: -14.4503     -74.1318        10.25
%         0    43    44    46    48    50    56
%         t-val = 4.1978 AAL: 48, coordinate: 10.4429     -84.4072     -13.6755
%         0    81    85
%         t-val = 3.2933 AAL: 85, coordinate: -51.6114      -42.358    -0.877889
%         86    90
%         t-val = 3.1861 AAL: 86, coordinate: 61.4449     -48.3575     -6.07904

%NEW
%    1    43    45    47    49    51    55
% t-val = 3.8645 AAL: 43, coordinate: -14.4503     -74.1318        10.25
%    1    81    85
% t-val = 3.4089 AAL: 81, coordinate: -64.7227     -25.7469      6.71894
%      1    43    44    46    48    50    56
% t-val = 4.1219 AAL: 48, coordinate: 10.4429     -84.4072     -13.6755
%      1    82    86    90
% t-val = 3.3386 AAL: 86, coordinate: 64.9292     -46.4401     -5.52821
    end
    
    
    
    %% presence versus perspective
    contrast = ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)))) - ... 
    ((((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) + ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)))) 
    slm  	   = SurfStatT(slm,contrast);
    statsprp   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.presence-perspective.FAC.days.'],1);
   
    
    numclus = sum(statsprp.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statsprp.clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statsprp.clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val 4.7795
% Cohen d 0.36103
% cluster sig 1.2527e-06
    end
    
    numclus = sum(statsprp.clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statsprp.clusid==i);
        unique(AAL2C2)
%         f= figure,
%         hist(AAL(statsprp.clusid==i),unique(AAL(statsprp.clusid==i)));
%         exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
        vertexid = find(statsprp.t == max(statsprp.t(statsprp.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
         maxtval = max(statsprp.t(statsprp.clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         0     4     6     8    10    14    16    20    24    26    34
%         t-val = 4.4333 AAL: 4, coordinate: 20.9469      21.1914      45.1848
%NEW
%      1     2     4     6     8    10    14    16    20    24    26    34
% t-val = 4.5074 AAL: 4, coordinate: 20.9469      21.1914      45.1848
    end

    numclus = sum(statsprp.inv_clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statsprp.inv_clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statsprp.inv_clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val -3.2869
% Cohen d -0.24829
% cluster sig 0.0050107
    end
    
    numclus = sum(statsprp.inv_clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statsprp.inv_clusid==i);
        unique(AAL2C2)
%         f= figure,
%         hist(AAL(statsprp.inv_clusid==i),unique(AAL(statsprp.inv_clusid==i)));
%         exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
        vertexid = find(statsprp.inv_t == max(statsprp.inv_t(statsprp.inv_clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statsprp.inv_t(statsprp.inv_clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         0    43    45    47    49    51
%         t-val = 3.2502 AAL: 47, coordinate: -8.21578      -85.058     -7.84128

%NEW
%      1    43    45    47    49    51
% t-val = 3.2648 AAL: 47, coordinate: -8.21578      -85.058     -7.84128
    end
    
    
    %% presence versus affect
    contrast   = ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)))) - ... 
    ((((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)) + ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1))))
    slm  	   = SurfStatT(slm,contrast);
    statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.presence-affect.FAC.days.'],1);
    
    numclus = sum(statspra.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statspra.clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statspra.clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val 4.84
% Cohen d 0.36561
% cluster sig 1.512e-06
% *******************************************************
% Post-hoc Clus: 2
% t-val 4.2122
% Cohen d 0.31818
% cluster sig 0.0016375
    end
    
    numclus = sum(statspra.clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statspra.clusid==i);
        unique(AAL2C2)
        f= figure,
        hist(AAL(statspra.clusid==i),unique(AAL(statspra.clusid==i)));
        exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
        vertexid = find(statspra.t == max(statspra.t(statspra.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statspra.t(statspra.clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         0     3     5     7     9    13    15    23    25    27
%         t-val = 3.5905 AAL: 3, coordinate: -32.4077      53.1346    -0.876738
%         0     6    10    16    26    28
%         t-val = 3.8899 AAL: 28, coordinate: 3.44955      43.2282     -21.7404

%NEW
%      1     3     5     7     9    13    15    23    25    27
% t-val = 3.5257 AAL: 3, coordinate: -32.4077      53.1346    -0.876738
%      1     6    10    16    26    28
% t-val = 3.7946 AAL: 28, coordinate: 3.44955      43.2282     -21.7404
    end

    numclus = sum(statspra.inv_clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,statspra.inv_clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(statspra.inv_clus.P(i))])
% *******************************************************
% Post-hoc Clus: 1
% t-val -5.5784
% Cohen d -0.42139
% cluster sig 1.2528e-06
% *******************************************************
% Post-hoc Clus: 2
% t-val -3.8816
% Cohen d -0.29322
% cluster sig 0.011749
% *******************************************************
% Post-hoc Clus: 3
% t-val -4.0282
% Cohen d -0.30428
% cluster sig 0.013263
    end
    
    numclus = sum(statspra.inv_clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(statspra.inv_clusid==i);
        unique(AAL2C2)
        f= figure,
        hist(AAL(statspra.inv_clusid==i),unique(AAL(statspra.inv_clusid==i)));
        exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
        vertexid = find(statspra.inv_t == max(statspra.inv_t(statspra.inv_clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(statspra.inv_t(statspra.inv_clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%         2     8    12    14    18    58
%         t-val = 3.9407 AAL: 2, coordinate: 44.9294       1.0018      46.7104
%         64    82
%         t-val = 3.4047 AAL: 82, coordinate: 64.0045     -33.1973      16.0255
%         0    18    58    64
%         t-val = 3.4129 AAL: 64, coordinate: 59.3674     -23.2183       31.263

%NEW
%      1     2     8    12    14    18    30    58    64    82
% t-val = 4.135 AAL: 2, coordinate: 47.1103      2.52934      46.6577
%     40    48    56
% t-val = 3.4587 AAL: 40, coordinate: 22.242     -23.3646     -26.1356
%      1    33    35    67    69
% t-val = 3.5479 AAL: 33, coordinate: -7.97467     -40.1174      38.1422
    end
    
    
    
    %% affect only
    contrast = ((((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)) + ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)))) - ...
         ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))))-...
        ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1))
    slm  	= SurfStatT(slm,contrast);
    statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.affect.FAC.days.'],1);
    
    %% affect G2
    contrast = ((((G4.Group_2.* DT.T3)-(G4.Group_2.* DT.T2)))) - ((((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)))) - ...
    ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0))))
    slm  	= SurfStatT(slm,contrast);
    statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.affectG2.FAC.days.'],1);
   
    %% affect G2
    contrast = ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0)) 
    slm  	= SurfStatT(slm,contrast);
    statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.affectG1.FAC.days.'],1);
   
   
    %% perspective only
    contrast = ((((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) + ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2))))- ...
        ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))))-...
        ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1))
    slm  	= SurfStatT(slm,contrast);
    statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.perspective.FAC.days.'],1);
    
    contrast = ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) - ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0))
    slm  	= SurfStatT(slm,contrast);
    statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.perspectiveG2.FAC.days.'],1);
   
     contrast = ((G4.Group_1.* DT.T3)-(G4.Group_1.* DT.T2)) - ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ...
         ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))
    slm  	= SurfStatT(slm,contrast);
    statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.perspectiveG1.FAC.days.'],1);
    
    %% presence only
    contrast = ((((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0)) + ((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))))
    slm  	= SurfStatT(slm,contrast);
    statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.presence.FAC.days.'],1);
   
    contrast = ((G4.Group_2.* DT.T1)-(G4.Group_2.* DT.T0))
    slm  	= SurfStatT(slm,contrast);
    statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.presenceG2.FAC.days.'],1);
    
    contrast = ((((G4.Group_1.* DT.T1)-(G4.Group_1.* DT.T0))))
    slm  	= SurfStatT(slm,contrast);
    statspra   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.025, ...
        [RPATH 'FIGURE2A.CT.presenceG1.FAC.days.'],1);
    
    %% post-hoc for significant seeds;
    numclus = sum(F5E_stats.inv_clus.P<0.05)
    for j = [1 : numclus]
    mysites={'Berlin','Leipzig'};
    for s = 1:2
        disp(mysites(s))
        k   = find(strcmp(site,mysites(s)));
        k   = intersect(keep,k);
        seedx       = mean(CT(k,find(F2E_stats.clusid==j)),2);
        seedxt      = term(seedx);
        ik1      = id(k,:);
        tfack1   = tp(k);
        g4k1     = group4(k);
        ak1      = age(k);
        sk1      = sex(k);
        
        A1       = term(ak1);
        S1      = term(sk1);
        G41      = term(g4k1);
        Sub1     = term(ik1);
        DT1      = term(tfack1);

        % model2
        M1       = 1 + Sub1 + A1 + S1 + G41 + DT1 + G41*DT1;
        slm = SurfStatLinMod(seedx,M1);
        slm = SurfStatT(slm, ((G41.Group_1.* DT1.T2)-(G41.Group_1.* DT1.T1)) - ((G41.Group_2.* DT1.T2)-(G41.Group_2.* DT1.T1)));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm.t)])
        slm = SurfStatT(slm, ((G41.Group_1.* DT1.T2)-(G41.Group_1.* DT1.T1))-((G41.Control.* DT1.T2)-(G41.Control.* DT1.T1)));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm.t)])
    end
    end
    
    %% sequence effect: what happens to previous significant seed changes? 
    numclus = sum(F1A_stats.clus.P<0.05)
    for j = [1 : numclus]
        seedx       = mean(CT(keep,find(F1A_stats.clusid==j)),2);
        seedxt      = term(seedx);
        slm = SurfStatLinMod(seedx,M);
        slm = SurfStatT(slm, ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm.t)])
        slm = SurfStatT(slm, ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)) - ((G4.Control.* DT.T2)-(G4.Control.* DT.T1)));
        disp([ 'TC1-RCC ROI: ' num2str(j) ',t= ' num2str(slm.t)])
        slm = SurfStatT(slm, ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)) - ((G4.Control.* DT.T2)-(G4.Control.* DT.T1)));
        disp([ 'TC2-RCC ROI: ' num2str(j) ',t= ' num2str(slm.t)])
        slm = SurfStatT(slm, ((G4.Group_1.* DT.T2)-(G4.Group_1.* DT.T1)));
        disp([ 'TC1 ROI: ' num2str(j) ',t= ' num2str(slm.t)])
        slm = SurfStatT(slm, ((G4.Group_2.* DT.T2)-(G4.Group_2.* DT.T1)));
        disp([ 'TC2 ROI: ' num2str(j) ',t= ' num2str(slm.t)])
        
%         ROI: 1,t= 0.056489
%         TC1-RCC ROI: 1,t= 1.5054
%         TC2-RCC ROI: 1,t= 1.4541
%         TC1 ROI: 1,t= 1.1438
%         TC2 ROI: 1,t= 1.0539
    end
    
    %% covariance affect > perspective
    j = [1]
    vertexid    = find(F2E_stats.t == max(F2E_stats.t(F2E_stats.clusid==j)))
    seedx       = mean(Ck(:,vertexid),2);
    seedxt      = term(seedx);
        %% covariance
    M       = 1 + Sub +  A + S + G4 + DT + G4*DT + G4*seedxt + DT*seedxt + DT*G4*seedxt;

    contrast = ((G4.Group_1.*DT.T2.*seedx)-(G4.Group_1.*DT.T1.*seedx))-((G4.Group_2.*DT.T2.*seedx)-(G4.Group_2.*DT.T1.*seedx))
    % FIGURE1A -- TC1+TC2
    slm = SurfStatLinMod(Ck,M,SW);
    slm  	= SurfStatT(slm,contrast);
    stats   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.05, ...
        [RPATH 'FIGURE1AA.CT.T2-T1.G1-G2.COV.peak.'],1);
    
    
    %% covariance affect > perspective
    j = [2]
    vertexid    = find(F2E_stats.inv_t == max(F2E_stats.inv_t(F2E_stats.inv_clusid==j)))
    seedx       = mean(Ck(:,vertexid),2);
    seedxt      = term(seedx);
        %% covariance
    M       = 1 + Sub +  A + S + G4 + DT + G4*DT + G4*seedxt + DT*seedxt + DT*G4*seedxt;

    contrast = ((G4.Group_2.*DT.T2.*seedx)-(G4.Group_2.*DT.T1.*seedx))-((G4.Group_1.*DT.T2.*seedx)-(G4.Group_1.*DT.T1.*seedx))
    % FIGURE1A -- TC1+TC2
    slm = SurfStatLinMod(Ck,M,SW);
    slm  	= SurfStatT(slm,contrast);
    stats   = BoSurfStatPlotAllStatsNew(slm, SM, mask, p_thresh, 0.05, ...
        [RPATH 'FIGURE1AA.CT.T2-T1.G2-G1.COV.peak.'],1);
    
end

%% ATTENTION X PRESENCE X TC3
for FIG3A = 1
    
    
    keep4   = find((strcmp(group4,'Group3')));
    keep2   = intersect(find(X(:,1) >-100), find(sum(X,2)~=0));
    size(mintersect(keep2,keep4, find(strcmp(groupN,'Affect')))) % 
    size(mintersect(keep2,keep4, find(strcmp(groupN,'Perspective')))) 
    size(mintersect(keep2,keep4, find(strcmp(groupN,'Presence')))) % N = 133
    keep3   = find(abs(mean(X,2))<0.1);
    size(mintersect(keep2,keep3, keep4)) % N = 132

    keep5   = find(abs(Att)<666);
    keep6   = find(attention>-666)
    keep7   = find(abs(Att_c)<666);
    keep8   = find(attention_c>-666);
    keep    = mintersect(keep2,keep3,keep4,keep5);

    size(keep) % N = 102
    size(mintersect(keep, find(strcmp(groupN,'Affect')))) 
    size(mintersect(keep, find(strcmp(groupN,'Perspective')))) 
    size(mintersect(keep, find(strcmp(groupN,'Presence')))) % N = 102
    mean(Att(keep)) % 0.05
    std(Att(keep))  % 0.15
    mask    = mean(CT,1)>0.4;
    Ck      = X(keep,:);
    ik      = id(keep,:);
    ink     = idnum(keep);
    tk      = tpnum(keep,:);
    tnk     = tp(keep,:);
    dk      = Days_last(keep);
    gk      = group(keep);
    g4k     = group4(keep);
    gNk     = groupN(keep);
    ak      = age(keep);
    sk      = sex(keep);
    Attk    = Att(keep);
    Attblk    = chAbl(keep);
    Att_ck    = Att_c(keep);
    Attbl_ck  = attention_c(keep);
    A       = term(ak);
    S       = term(sk);
    G       = term(gk);
    G4      = term(g4k);
    GN      = term(gNk);
    D       = term(dk);
    Sub     = term(ik);
    T       = term(tk);
    Tn      = term(tnk);
    ATT     = term(Attk);
    ATTBL   = term(Attblk);
    ATTc     = term(Att_ck);
    ATTBLc   = term(Attbl_ck);
    CTM     = mean(Ck,2);
    CH      = term(CTM);

    M       = 1  + GN + A + S + random(Sub) + I;
    slm = SurfStatLinMod(Attk,M);
    slm = SurfStatT(slm, (GN.Presence)-(GN.Affect)-(GN.Perspective));
    disp([ ',t= ' num2str(slm.t)])
    % 2.48

    %% whole brain attention change; 
    M       = 1 + CH + A + S  + GN + ATT + GN*ATT + random(Sub) + I;
    slm     = SurfStatLinModS(Ck,M,SW);
    slm  	= SurfStatT(slm,(GN.Affect.*Attk));
    
    f=figure,
    BoSurfStatViewData(slm.t,SInf,'')
    SurfStatColLim([1.69 5])
    colormap([0 0 0; mycol.blackblue])
    exportfigbo(f,[RPATH 'ATT_mod.png'],'png',10)
    clusatt    = BoSurfStatPlotAllStats(slm,SM, mask, 0.025, 0.025, [RPATH '/FIG1C.presence.T3.att_CHk.WB' ],1);
    clusatt    = BoSurfStatPlotAllStats(slm,SInf, mask, 0.025, 0.025, [RPATH '/FIG1C.presence.T3.inf.att_CHk.WB' ],1);

    numclus = sum(clusatt.clus.P < 0.025);
    for i = 1 : numclus
        seedx       = mean(CT(keep,clusatt.clusid==i),2);
        seedxt      = term(seedx);
        
        M2       =1 + CH + A + S + GN + ATT + GN*ATT + Sub;
        slm1 = SurfStatLinMod(seedx,M2);
        slm1 = SurfStatT(slm1, (GN.Presence.*Attk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
% ROI: 0+1i,t= 2.6306
        M3       = CH + A + S + (1 + ATT)*(1 + GN)*(1 + random(Sub)) + I;
        slm1 = SurfStatLinMod(seedx,M3);
        slm1 = SurfStatT(slm1, (GN.Presence.*Attk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
% ROI: 0+1i,t= 2.0724
    end
    
    contrast = (GN.Presence.*Attk)
    numclus = sum(clusatt.clus.P<0.025)
    for i = 1:numclus
        ROI = mean(CT(keep,clusatt.clusid==i),2);
        slma      = SurfStatLinMod(ROI,M);
        slma      = SurfStatT(slma,contrast); 
        disp('*******************************************************')
        disp(['Post-hoc Clus: ' num2str(i)]) 
        disp([ 't-val ' num2str(slma.t)])   
        d = 2*(slma.t) / sqrt(slma.df);
        disp([ 'Cohen d ' num2str(d)])
        disp([ 'cluster sig ' num2str(clusatt.clus.P(i))])
%         *******************************************************
%         Post-hoc Clus: 1
%         t-val 2.6306
%         Cohen d 0.30949
%         cluster sig 2.2688e-06
    end
    
    numclus = sum(clusatt.clus.P < 0.025);
    for i = 1 : numclus
        AAL2C2 = AAL(clusatt.clusid==i);
        unique(AAL2C2)
        f= figure,
        hist(AAL(clusatt.clusid==i),unique(AAL(clusatt.clusid==i)));
        exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
        vertexid = find(clusatt.t == max(clusatt.t(clusatt.clusid==i)));
        aalnum = AAL(:,vertexid);
        coordinate = SM.coord(:,vertexid)';
        maxtval = max(clusatt.t(clusatt.clusid==i));
        disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%     0    51    55    65    79    81    85    89
%     t-val = 4.205 AAL: 85, coordinate: -64.5411     -35.2448     -14.0829

% NEW
%      1    55    81    85    89
% t-val = 3.5417 AAL: 89, coordinate: -58.2058     -32.3544     -18.8977
    end
    
    
    %% attention
    numclus = sum(clusatt.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(clusatt.clusid==j)),2);
         M    = 1 + CH + A + S; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepatt1  = find(abs(Att)<666);
         keepp = find(strcmp(groupN,'Presence'));
         keepatt = mintersect(keep,keepatt1,keepp);
         f= figure,
         hold on
         scatter(Att(keepatt),targetresid(keepatt), 60, 'Marker','o', 'MarkerFaceColor', [0.9 0.9 0], 'MarkerEdgeColor', [0 0 0]), lsline
         %axis([-0.3 0.3 -1 1])
         exportfigbo(f,[RPATH 'FIG3C.Attch.Presence.scatter' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepatt),Att(keepatt))
%          r = 0.34
    end
    
    numclus = sum(clusatt.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(clusatt.clusid==j)),2);
         M    = 1; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepc1  = find(abs(Comp)<666);
         keepp = find(strcmp(groupN,'Presence'));
         keepc = mintersect(keep,keepc1,keepp);
         f= figure,
         hold on
         scatter(Comp(keepc),targetresid(keepc), 60, 'Marker','o', 'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0 0 0]), lsline
         exportfigbo(f,[RPATH 'FIG3C.Comp.PresenceCH.scatter' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepc),Comp(keepc))
         % r = -0.04   p = 0.73
    end
    
    for j = [1 : numclus]
         target = mean(Ck(:,find(clusatt.clusid==j)),2);
         M    = 1 ; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepc1  = find(abs(Tom)<666);
         keepp = find(strcmp(groupN,'Presence'));
         keepp = mintersect(keep,keepc1,keepp);
         f= figure,
         hold on
         scatter(Tom(keepp),targetresid(keepp), 60, 'Marker','o', 'MarkerFaceColor', [0 0.7 0],'MarkerEdgeColor', [0 0 0]), lsline
         %axis([-0.3 0.3 -1 1])
         exportfigbo(f,[RPATH 'FIG3C.ToM.PresenceCh.scatter' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepp),Tom(keepp))
%          r = -0.0913 p =0.36
    end
  
    
    %% separate groups 
    
    numclus = sum(clusatt.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(clusatt.clusid==j)),2);
         M    = 1 + CH + A + S; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepatt1  = find(abs(Att)<666);
         keepp = find(strcmp(groupN,'Presence'));
         keepp2 = find(strcmp(group,'Group_1'));
         keepatt = mintersect(keep,keepatt1,keepp,keepp2);
         f= figure,
         hold on
         scatter(Att(keepatt),targetresid(keepatt), 60, 'Marker','o', 'MarkerFaceColor', [0.9 0.9 0], 'MarkerEdgeColor', [0 0 0]), lsline
         axis([-0.3 0.6 -0.2 0.2])
         exportfigbo(f,[RPATH 'FIG3C.Attch.Presence.scatter.G1.' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepatt),Att(keepatt))
%          r = 0.45   p =   0.0007
    end
    
    numclus = sum(clusatt.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(clusatt.clusid==j)),2);
         M    = 1 + CH + A + S; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepatt1  = find(abs(Att)<666);
         keepp = find(strcmp(groupN,'Presence'));
         keepp2 = find(strcmp(group,'Group_2'));
         keepatt = mintersect(keep,keepatt1,keepp,keepp2);
         f= figure,
         hold on
         scatter(Att(keepatt),targetresid(keepatt), 60, 'Marker','o', 'MarkerFaceColor', [0.9 0.9 0], 'MarkerEdgeColor', [0 0 0]), lsline
         axis([-0.3 0.6 -0.2 0.2])
         exportfigbo(f,[RPATH 'FIG3C.Attch.Presence.scatter.G2.' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepatt),Att(keepatt))
%          r = 0.17
    end
    
    
        
    numclus = sum(F1A_stats.peak.P<0.025)
    for j = [1]
%         seedx       = mean(Ck(:,F1A_stats.peak.vertid(j)),2);
%         seedxt      = term(seedx);
        mymask = zeros(1,size(CT,2));
        mymask(19376) = 1;
        mymask = SurfStatSmooth(mymask,SW,10);
        mymask = mymask / max(mymask);
        mymask = mymask > 0.5;
        seedx       = mean(Ck(:,mymask),2);
        seedxt      = term(seedx);
        f = figure,
        BoSurfStatViewData(mymask,SM,'')
        exportfigbo(f,[RPATH 'FIG.presence.png'],'png',10)
       %% simple ROI
        M       = 1  +  GN + A + S + ATT + GN*ATT + random(Sub) + I;
        slm = SurfStatLinModS(seedx,M);
        slm = SurfStatT(slm, (GN.Presence.*Attk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm.t)])
        % ROI: 1,t= -0.43645
        M       = CH + A + S + (1 + ATT)*(1 + GN)*(1 + random(Sub)) + I;
        slm = SurfStatLinModS(seedx,M);
        slm = SurfStatT(slm, (GN.Presence.*Attk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm.t)])
        % ROI: 1,t= -0.97542
        %% covariance
        M       = 1  + GN + A + S + seedxt+ GN*seedxt + CH + random(Sub) + I;
        slm         = SurfStatLinModS(Ck,M,SW);
        slm         = SurfStatT(slm, (GN.Presence.*seedx));
        clus1C_1    = BoSurfStatPlotAllStatsNew(slm,SM, mask, 0.025, 0.025, [RPATH '/FIG1C.presnece.cov.clusid.APRIORI.AI' ],1);
        

       %% modulation
        M       = 1 + CH + GN + A + S + ATT + seedxt + GN*ATT + GN*seedxt + ATT*seedxt + ATT*GN*seedxt + random(Sub) + I;
        slm         = SurfStatLinModS(Ck,M,SW);
        slm         = SurfStatT(slm, (GN.Presence.*Attk.*seedx));
        clus1C_2    = BoSurfStatPlotAllStats(slm,SM, mask, 0.025, 0.025, [RPATH '/FIG1C.presence.cov.att.clusid' ],1);
    
        contrast = (GN.Presence.*Attk.*seedx)
        numclus = sum(clus1C_2.clus.P<0.025)
        for i = 1:numclus
            ROI = mean(CT(keep,clus1C_2.clusid==i),2);
            slma      = SurfStatLinMod(ROI,M);
            slma      = SurfStatT(slma,contrast);
            disp('*******************************************************')
            disp(['Post-hoc Clus: ' num2str(i)])
            disp([ 't-val ' num2str(slma.t)])
            d = 2*(slma.t) / sqrt(slma.df);
            disp([ 'Cohen d ' num2str(d)])
            disp([ 'cluster sig ' num2str(clus1C_2.clus.P(i))])
%             Post-hoc Clus: 1
%             t-val 2.994
%             Cohen d 0.36714
%             cluster sig 0.0024627
        end
        
        numclus = sum(clus1C_2.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(clus1C_2.clusid==i);
            unique(AAL2C2)
            f= figure,
            hist(AAL(clus1C_2.clusid==i),unique(AAL(clus1C_2.clusid==i)));
            exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
            vertexid = find(clus1C_2.t == max(clus1C_2.t(clus1C_2.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            disp([ 'AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%    1     4     6     8    10    24    26
%    AAL: 8, coordinate: 35.3817      52.0318        7.687
        end
        
    
    end
    
        seedx       = mean(Ck(:,inter_str==1),2);
        seedxt      = term(seedx);
    
        %% simple ROI
        M       = 1  + CH + GN + A + S + ATT + GN*ATT + random(Sub) + I;
        slm = SurfStatLinMod(seedx,M);
        slm = SurfStatT(slm, (GN.Presence.*Attk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm.t)])
        % ROI: 1,t= -1.1703
    
         target = Attk
         M    = 1 + ATTBL; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         Attcorrk = targetresid(keep,:)
         Attck  = term(Attcorrk)
             
  
    
end

%% COMPAsSSION x AFFECT X TC3
for FIG3B = 1
    

    keep1   = find((strcmp(group4,'Group3')));
    keep2   = intersect(find(X(:,1) >-666), find(sum(X,2)~=0));
    keep3   = find(abs(mean(X,2))<0.1);
    keep4   = find((strcmp(group4,'Group3')));
    size(mintersect(keep2,keep3, keep4)) % N = 126 (127)
    % odd insula
    keep5   = find((strcmp(id,'22085')));
    mask    = mean(CT,1)>0.9;
    keep7   = find(compassion>0);
    keep7b  = find(abs(Comp)<666);
    keep    = mintersect(keep1, keep2,keep7,keep7b, keep3);
    
    size(keep) 
    size(mintersect(keep, find(strcmp(groupN,'Affect')))) % N = 121
    size(mintersect(keep, find(strcmp(groupN,'Perspective')))) 
    size(mintersect(keep, find(strcmp(groupN,'Presence')))) 

    Ck      = X(keep,:);
    ik      = id(keep,:);
    ink     = idnum(keep);
    tk      = tpnum(keep,:);
    tnk     = tp(keep,:);
    gk      = group(keep);
    g4k     = group4(keep);
    gNk     = groupN(keep);
    ak      = age(keep);
    sk      = sex(keep);
    Compk   = Comp(keep);
    dk      = Days_last(keep);
    combf   = chcbl(keep);
    A       = term(ak);
    S       = term(sk);
    G       = term(gk);
    G4      = term(g4k);
    GN      = term(gNk);
    Sub     = term(ik);
    CMP     = term(Compk);
    T       = term(tk);
    Tn      = term(tnk);
    D       = term(dk);
    
    CHk   = mean(Ck,2);
    CH    = term(CHk);
    
    Cbf   = term(combf);
    
    Tbf   = mean(CT(keep,:),2);
    TF    = term(Tbf);
    %% correlation before and change
    M       = 1  + GN + A + S + random(Sub) + I;
    slm = SurfStatLinMod(Compk,M);
    slm = SurfStatT(slm, (GN.Affect)-(GN.Presence)-(GN.Perspective));
    disp([ ',t= ' num2str(slm.t)])
    % 3.36
    M       = 1  + GN + A + S + random(Sub) + I;
    slm = SurfStatLinMod(Compk,M);
    slm = SurfStatT(slm, (GN.Presence));
    disp([ ',t= ' num2str(slm.t)])
    % -3.65
    M       = 1  + GN + A + S + random(Sub) + I;
    slm = SurfStatLinMod(Compk,M);
    slm = SurfStatT(slm, (GN.Perspective));
    disp([ ',t= ' num2str(slm.t)])
    % 0.76 
    
    target = Compk
         M    = 1 + Cbf; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         Comcorrk = targetresid(keep,:)
         Comck  = term(Comcorrk)
         
    var = Comp(keep);
    var_fac = sk;
    var_fac(find(var<= median(var))) = {'low'};
    var_fac(find(var> median(var))) = {'high'};
    vf = term(var_fac);
    
    %% whole brain
    M       = 1 + CH + A + S + GN + CMP + GN*CMP + random(Sub) + I;
    slm     = SurfStatLinModS(Ck,M,SW);
    slm  	= SurfStatT(slm,(GN.Affect.*Compk)-(GN.Presence.*Compk));
    
    f=figure,
    BoSurfStatViewData(slm.t,SInf,'')
    SurfStatColLim([1.69 5])
    colormap([0 0 0; mycol.blackblue])
    exportfigbo(f,[RPATH 'FIG3C.CMP_blank.png'],'png',10)
    cluscomp   = BoSurfStatPlotAllStats(slm,SM, mask, 0.025, 0.025, [RPATH '/FIG3C.affect-presence.TC3.comp.WB' ],1);
    cluscomp   = BoSurfStatPlotAllStats(slm,SInf, mask, 0.025, 0.025, [RPATH '/FIG3C.affect.INF.comp.WB' ],1);

   %% test for complex model with random effects
    numclus = sum(cluscomp.clus.P < 0.025);
    for i = 1 : numclus
        seedx       = mean(CT(keep,cluscomp.clusid==i),2);
        seedxt      = term(seedx);
        
        M2       = 1 + CH + A + S + GN + CMP + GN*CMP + Sub;
        slm1 = SurfStatLinMod(seedx,M2);
        slm1 = SurfStatT(slm1, (GN.Affect.*Compk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
        
        M3       = CH + A + S + (1 + CMP)*(1 + GN)*(1 + random(Sub)) + I;
        slm1 = SurfStatLinMod(seedx,M3);
        slm1 = SurfStatT(slm1, (GN.Affect.*Compk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm1.t)])
        
% ROI: 1,t= 3.05
% ROI: 1,t= 2.04
    end
    
    
    
    
        contrast = (GN.Affect.*Compk)
        numclus = sum(cluscomp.clus.P<0.025)
        for i = 1:numclus
            ROI = mean(CT(keep,cluscomp.clusid==i),2);
            slma      = SurfStatLinMod(ROI,M);
            slma      = SurfStatT(slma,contrast);
            disp('*******************************************************')
            disp(['Post-hoc Clus: ' num2str(i)])
            disp([ 't-val ' num2str(slma.t)])
            d = 2*(slma.t) / sqrt(slma.df);
            disp([ 'Cohen d ' num2str(d)])
            disp([ 'cluster sig ' num2str(cluscomp.clus.P(i))])
            % *******************************************************
            % Post-hoc Clus: 1
            % t-val 2.9703
            % Cohen d 0.31529
            % cluster sig 2.0786e-06
        end
        
        numclus = sum(cluscomp.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(cluscomp.clusid==i);
            unique(AAL2C2)
            f= figure,
            hist(AAL(cluscomp.clusid==i),unique(AAL(cluscomp.clusid==i)));
            exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
            vertexid = find(cluscomp.t == max(cluscomp.t(cluscomp.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            maxtval = max(cluscomp.t(cluscomp.clusid==i));
            disp([ 't-val = ' num2str(maxtval) ' AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
%           0    12    14    16    18    30    82    84
%           t-val = 3.6604 AAL: 30, coordinate: 31.9174       18.878      7.75285

%NEW
%     1    12    14    16    30
% t-val = 3.9006 AAL: 14, coordinate: 32.1758      22.3563      9.10549

        end
        
    %% attention
    numclus = sum(cluscomp.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(cluscomp.clusid==j)),2);
         M    = 1 +  CH + A + S; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepatt1  = find(abs(Att)<666);
         keepp = find(strcmp(groupN,'Affect'));
         keepatt = mintersect(keep,keepatt1,keepp);
         f= figure,
         hold on
         scatter(Att(keepatt),targetresid(keepatt),60, 'Marker','o', 'MarkerFaceColor', [0.9 0.9 0], 'MarkerEdgeColor', [0 0 0.2]), lsline
         %axis([-0.5 0.3 -0.4 0.3])
         exportfigbo(f,[RPATH 'FIG3C.Attch.affect.scatter' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepatt),Att(keepatt))
         % r = -0.09 , p<0.42
    end
    
    numclus = sum(cluscomp.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(cluscomp.clusid==j)),2);
         M    = 1 + CH; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepc1  = find(abs(Comp)<666);
         keepp = find(strcmp(groupN,'Affect'));
         keepc = mintersect(keep,keepc1,keepp);
         f= figure,
         hold on
         scatter(Comp(keepc),targetresid(keepc),60,'Marker','o', 'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0 0 0]), lsline
         %axis([-1.5 2.5 -0.65 0.4])
         exportfigbo(f,[RPATH 'FIG3C.Comp.affect.scatter' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepc),Comp(keepc))
         %r =0.23
    end
    
    for j = [1 : numclus]
         target = mean(Ck(:,find(cluscomp.clusid==j)),2)
         M    = 1 ; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepc1  = find(abs(Tom)<666);
         keepp = find(strcmp(groupN,'Affect'));
         keepp = mintersect(keep,keepc1,keepp);
         f= figure,
         hold on
         scatter(Tom(keepp),targetresid(keepp), 60, 'Marker','o', 'MarkerFaceColor', [0 0.7 0], 'MarkerEdgeColor', [0 0 0]), lsline
         %axis([-0.3 0.3 -1 1])
         exportfigbo(f,[RPATH 'FIG3C.ToM.affect.scatter' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepp),Tom(keepp))
         % r = 0.0934, p < 0.3082
    end
    
    
    %% separate groups
    numclus = sum(cluscomp.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(cluscomp.clusid==j)),2);
         M    = 1; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepc1  = find(abs(Comp)<666);
         keepp = find(strcmp(groupN,'Affect'));
         keepp2 = find(strcmp(groupN,'Group3'));
         keepc = mintersect(keep,keepc1,keepp2);
         f= figure,
         hold on
         scatter(Comp(keepc),targetresid(keepc),60,'Marker','o', 'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0 0 0]), lsline
         axis([-1.5 2.5 -0.2 0.2])
         exportfigbo(f,[RPATH 'FIG3C.Comp.affect.scatter.G1.' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepc),Comp(keepc))
         %r =0.31, p<0.02
    end
    
    numclus = sum(cluscomp.clus.P<0.025)
    for j = [1 : numclus]
         target = mean(Ck(:,find(cluscomp.clusid==j)),2);
         M    = 1 + CH + A + S ; 
         slm = SurfStatLinMod(target,M); 
         targetresid(keep,:) = target - slm.X*slm.coef; 
         keepc1  = find(abs(Comp)<666);
         keepp = find(strcmp(groupN,'Affect'));
         keepp2 = find(strcmp(group,'Group_2'));
         keepc = mintersect(keep,keepc1,keepp,keepp2);
         f= figure,
         hold on
         scatter(Comp(keepc),targetresid(keepc),60,'Marker','o', 'MarkerFaceColor', [0.8 0 0], 'MarkerEdgeColor', [0 0 0]), lsline
         axis([-1.5 2.5 -0.2 0.2])
         exportfigbo(f,[RPATH 'FIG3C.Comp.affect.scatter.G2.' num2str(j) '.png'],'png',10)
         [r p] = corr(targetresid(keepc),Comp(keepc))
         %r =0.21
    end
    
    
    
    numclus = sum(F2E_stats.clus.P<0.025)
    for j = [1 : numclus]
        mymask = zeros(1,size(CT,2));
        mymask(13269) = 1;
        mymask = SurfStatSmooth(mymask,SW,10);
        mymask = mymask / max(mymask);
        mymask = mymask > 0.5;
        seedx       = mean(X(keep,mymask),2);
        seedxt      = term(seedx);
        
        f= figure,
        BoSurfStatViewData(mymask,SM,'')
        exportfigbo(f,[RPATH 'FIG.AI.png'],'png',10)

        %% simple ROI
        M       = 1  + CH + GN + A + S + CMP + GN*CMP + random(Sub) + I;
        slm = SurfStatLinMod(seedx,M);
        slm = SurfStatT(slm, (GN.Affect.*Compk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm.t)])
        % pCC = -0.97
        % dlPFC = 1.10
        % oper = 0.002
        % mCC = -1.3
        % AI = 0.87
        %% covariance
        M       = 1  + CH + GN + A + S + CMP + seedxt + GN*CMP + GN*seedxt + CMP*seedxt + CMP*GN*seedxt + random(Sub) + I;
        slm         = SurfStatLinModS(Ck,M,SW);
        slm         = SurfStatT(slm, (GN.Affect.*Compk.*seedx));
        clusempC    = BoSurfStatPlotAllStats(slm,SM, mask, 0.025, 0.025, [RPATH '/FIG3C.affect.cov.CMP.AI.'  ],1);
        
        contrast = (GN.Affect.*Compk.*seedx)
        numclus = sum(clusempC.clus.P<0.025)
        for i = 1:numclus
            ROI = mean(CT(keep,clusempC.clusid==i),2);
            slma      = SurfStatLinMod(ROI,M);
            slma      = SurfStatT(slma,contrast);
            disp('*******************************************************')
            disp(['Post-hoc Clus: ' num2str(i)])
            disp([ 't-val ' num2str(slma.t)])
            d = 2*(slma.t) / sqrt(slma.df);
            disp([ 'Cohen d ' num2str(d)])
            disp([ 'cluster sig ' num2str(clusempC.clus.P(i))])
%             *******************************************************
%             Post-hoc Clus: 2
%             t-val 3.1825
%             Cohen d 0.34932
%             cluster sig 0.016796
        end
        
        numclus = sum(clusempC.clus.P < 0.025);
        for i = 1 : numclus
            AAL2C2 = AAL(clusempC.clusid==i);
            unique(AAL2C2)
            f= figure,
            hist(AAL(clusempC.clusid==i),unique(AAL(clusempC.clusid==i)));
            exportfigbo(f,[RPATH 'G2presence' num2str(i) '.png'],'png',10);
            vertexid = find(clusempC.t == max(clusempC.t(clusempC.clusid==i)));
            aalnum = AAL(:,vertexid);
            coordinate = SM.coord(:,vertexid)';
            disp([ 'AAL: ' num2str(aalnum) ', coordinate: '  num2str(coordinate)])
            %   4     8    24
            %   AAL: 4, coordinate: 14.6581      43.2551      41.7543
        end
    end
    
    %% simple correlation
    numclus = sum(affectiveF2.stats.clus.P<0.025)
    for j = [1 : numclus]
        seedx       = mean(X(keep,affectiveF2.stats.clusid==j),2);
        seedxt      = term(seedx);
        
        f= figure,
        BoSurfStatViewData(affectiveF2.stats.clusid==2,SM,'')


        %% simple ROI
        M       = 1  + CH + GN + A + S + CMP + GN*CMP + Sub;
        slm = SurfStatLinMod(seedx,M);
        slm = SurfStatT(slm, (GN.Affect.*Compk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm.t)])
    end
    
    
    numclus = sum(T2T0stats.clus.P<0.025)
    for j = [1 : numclus]
        seedx       = mean(X(keep,T2T0stats.clusid==j),2);
        seedxt      = term(seedx);
        
        f= figure,
        BoSurfStatViewData(T2T0stats.clusid==j,SM,'')


        %% simple ROI
        M       = 1  + CH + GN + A + S + CMP + GN*CMP + Sub;
        slm = SurfStatLinMod(seedx,M);
        slm = SurfStatT(slm, (GN.Affect.*Compk));
        disp([ 'ROI: ' num2str(j) ',t= ' num2str(slm.t)])
    end
end








  

    
    
    
    

  
  
   
 
    
