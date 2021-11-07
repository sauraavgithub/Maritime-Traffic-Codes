% Data Extraction from the AIS data


clear all;
clc;

% 1
% Defining the area or Zone of analysis
% Enter the Longitude and Latitude (geographical coordinates)in arrays for the zone to be analysed.

lon = [69.0214735 72.744735 72.900597 69.370597]; % Longitude
lat = [19.489866 19.489866 18.490359 18.490359];           % Latitude
 
% 2
% Applying UTM (Universal Transverse Mercator) Coordinate transformation system

p1 = [lat(1),lon(1)];
z1 = utmzone(p1);
[ellipsoid,estr] = utmgeoid(z1);
utmstruct = defaultm('utm'); 
utmstruct.zone = z1; 
utmstruct.geoid = ellipsoid(1,:); 
utmstruct = defaultm(utmstruct); 
[utm_x utm_y]= mfwdtran(utmstruct,lat,lon);

%3
%Reading all *.mat files

f_data = dir('Mat_files/*.mat');
N = length(f_data); % Number of .mat files
number = 0;
     
load('V_GT.mat'); % Contains two dimensional array of MMSI & GT for various Vessels (vessel databse)  
count0 = length(V_GT);
     
     for f_count = 1:N
                 
        filename=['Mat_files\',f_data(f_count).name];
        %Data = importdata(filename);
        load(filename);
        
        l = length(f_data(f_count).name);
        t_indx = f_data(f_count).name(5:l-4); 
        if size(data_extract) == [0,0]
            for row = 1:361
                united_CR(row,1) = 0;
                all_CR(row,1) = 0;
            end
%             f_name_1 = ['Collision Risk\','zone\','Risk',t_indx,'.mat']
%             save(f_name_1,'united_CR');
%               f_name_1 = ['CR\','zone\','Risk',t_indx,'.mat']
%               save(f_name_1,'all_CR');
        else
        time = [];
        %Reading data from each .mat file into variables.
        %Date = Data(:,1); includes both colum A and B of excel file.
        %[Y,M,D,HH,MM,SS]= datevec(Date);
        HH = data_extract(:,2);
        MM = data_extract(:,3);
        SS = data_extract(:,4);
        ll = length(HH);
        for nn = 1:ll
            temp_t = HH(nn)*3600 + MM(nn)*60 + SS(nn);
            time = [time; temp_t];
        end
%         display(time);
        vessel_id = data_extract(:,5);
        SOG = data_extract(:,9);
        longitude = data_extract(:,7);
        Latitude = data_extract(:,6);
        COG = data_extract(:,10);
 
        u_vessel_id = unique(vessel_id); 
           
        count1 =length(u_vessel_id);
        
        %4
        %Filtering AIS data for particular zone based on geographical coordinates for analysis
        
        indx_lon = find(longitude>= min(lon) & longitude<= max(lon));
        indx_lat = find(Latitude>=min(lat) & Latitude<=max(lat));
        indx = intersect(indx_lon,indx_lat);
        
        % Extracting AIS data for each vessel corresponding to their MMSI no.
         
         for i=1:count1 
             
                msg = strcat('file_',num2str(f_count),'_',num2str(count1),'_',num2str(i));
                display(msg);
                
                a_chk = find(V_GT(:,1)== u_vessel_id(i,1)); %Cross checking with the vessel database
                display(a_chk);
                
            if (a_chk ~=0)
                
                MMSI = num2str(u_vessel_id(i,1));
                Data_extract=[];
            
                if  length(MMSI)==9  % MMSI no. of vessels are of 9 digits
                    
                    indx_mmsi = find(vessel_id== str2num(MMSI));
                    indx1 = intersect(indx_mmsi,indx);
                    count2 = length(indx1);
                    
                    for j = 1:count2
                        tf=0;
                        p = indx1(j,1);     
                        Lon = longitude(p,1);
                        Lat = Latitude(p,1);
                        [x y]= mfwdtran(utmstruct,Lat,Lon);
                
                        
                        % Checking whether (x,y) lies within the zone defined by [utm_x,utm_y]
                        
                        z_chk = zone_chk(x,y,utm_x,utm_y);
                   
            
                        if z_chk == 1 
                            
                            number = number + 1;
                    
                            if SOG(p,1)==102.3 %% vague value of speed when transmission is absent mostly when vessel is at rest.
                    
                                SOG(p,1)=0; 
                               
                            end
                            
                            
                            % Checking for Vessels permissible velocity for analysis
                            
                            v_chk = velocity_chk(SOG(p,1));
                 
                            if v_chk ==1 
            
                                A1 = [time(p,1) longitude(p,1) Latitude(p,1) COG(p,1) SOG(p,1)];
        
                                C=isnan(A1); %% checking for 'NaN' entries.
                                ind = find(C==1);
                                
                                if isempty(ind)
                                                           
                                    tf=1;                    
                                    
                                end
        
                                if tf==1
                                     
                                    temp = [time(p,1) Latitude(p,1) longitude(p,1) COG(p,1) SOG(p,1) ];
                                    
                                    % Storing usable AIS information in an array (Data_extract) corresponding to each MMSI no.
                                    if isempty(Data_extract)
                                        
                                       Data_extract = temp;
                                    else
                                        
                                       Data_extract = [Data_extract;temp]; 
                                    end
                                    
                                    % Storing the extracted AIS data for each vessel in a file named in the format 'MMSI'.mat
                                    
                                    fname = ['output\','zone\',MMSI,'.mat'];
                                    save(fname,'Data_extract');                               
                                    
                                    
                                    load(fname);
                                    A2 = Data_extract;
                                    
                                    if  size(A2,1)>2
                                        minTime = A2(1,1);
                                        maxTime = A2(end,1);
                                        xi= (minTime:maxTime)';
                                        
                                        % Interpolating the extracted AIS data for each 'MMSI'.mat file to generate a
                                        % continuous time series for analysis.
                                        Data_sync = interpolateMatrix(minTime,maxTime,[A2(:,1) A2(:,2:5)]);
    
                                        fname = ['output_sync\','zone\',MMSI,'.mat'];
                                        save(fname,'Data_sync'); %Data_sync implies interpolated Data_extract.
                                    end
                                end

                            end
                        end
                    end

                end
            end
         end
lon = [69.0214735 72.744735 72.900597 69.370597];
lat = [19.489866 19.489866 18.490359 18.490359];
         

p1 = [lat(1),lon(1)];
z1 = utmzone(p1);
[ellipsoid,estr] = utmgeoid(z1);
utmstruct = defaultm('utm'); 
utmstruct.zone = z1; 
utmstruct.geoid = ellipsoid(1,:); 
utmstruct = defaultm(utmstruct); 
[utm_x utm_y]= mfwdtran(utmstruct,lat,lon);

%Reading all filenames from folder output_sync.
f_name = dir(fullfile('output_sync/zone/*.mat'));
N = length(f_name);


% united_CR=[];
% all_CR=[];
interval = 10; % time interval in seconds at which collision risk value is calculated 
p=0;
Matrx = [];

% C is the visibility criteria(Day/Night). Depending upon the t_indx, value of C is decided.

t_indx1 = mod(str2double(t_indx),24);
if (t_indx1<=17 && t_indx1>=8)
    C=1; % Day condition
else
    C=2; % Nigt condition
end

if  t_indx1==0
    t_indx1=str2double(t_indx);
end

time = 3600*(t_indx1-1);

for T = time:time+3600 % T varies for a duration of 3600 secconds (1 hour).
    
    if rem(T,interval)==0 % Collision risk value is calculated at pre-defined interval.
       k = 0;
       fname1={};
       
       % Checking no. of ships in the zone present at time T.
       for count = 1:N
    
            fname = ['output_sync\zone\',f_name(count).name];
            load(fname);
    
            x1 = find(Data_sync(:,1)== T);   
    
            tf =isempty(x1);
            if tf~=1
                k=k+1;
                %fname1 stores the filename and MMSI no. for vessels
                %present in the zone at time T.
                 fname1(k,1:2) = {fname,char(f_name(count).name(1:9))}; 
            end
     
       end

        count1 = k; % No. of ships in the zone present at time T.

        CR = zeros(count1,count1);
        dist=zeros(count1,count1);
        
        if count1~=0
            
            %Carrying our pair wise interaction of ships to generate collision risk matrix [CR]
            
            for i = 1: count1 
                
                    load(char(fname1(i,1)));
                    
                    mmsi = str2double(fname1(i,2));
                    indx_GT = find(V_GT == mmsi);
                    VC  = V_GT(indx_GT,2);
                    
                    ind0 = find(Data_sync(:,1)== T);
                    pos0_x = Data_sync(ind0,2); pos0_y = Data_sync(ind0,3); c0 = deg2rad(Data_sync(ind0,4)); v0 = Data_sync(ind0,5)*0.5144;
                    [x0 y0] = mfwdtran(utmstruct,pos0_y,pos0_x);
                    
                    % Reading the coefficient values.
%                     [beta thresholds risk_score] =  coefficients(VC);
    
                    for j = 1:count1
        
                    if (i==j)
%                         CR(i,j)=0;
                        dist(i,j)=0;
                    else
                        load(char(fname1(j,1)));
                        ind1 = find(Data_sync(:,1)== T);
                        pos_x = Data_sync(ind1,2); pos_y = Data_sync(ind1,3); c = deg2rad(Data_sync(ind1,4)); v = Data_sync(ind1,5)*0.5144; 
                        [x y] = mfwdtran(utmstruct,pos_y,pos_x);
         
                        dist(i,j) = Dist(x0,y0,x,y); %Distance between two points.
                        display (dist);
                        b_rt = RB(x0,y0,x,y,c0); % Relative Bearing
                        v_rt = sqrt(v0^2+ v^2-2*v0*v*cos(c-c0)); % Relative Velocity
                        c_rt = acos((v0-v*cos(c-c0))/v_rt); % Relative course
                         
                        % Calculating the proximity parameters.
                        [DCPA,TCPA]= dcpa1(dist(i,j),c_rt,b_rt,v_rt);
                        matrix = [DCPA,TCPA];
                        Matrx = [Matrx;matrix];
                        

%                         if TCPA>=0
%                             Prob=zeros(1,5);
%                             for m=1:5 % for five risk criteria (VHR, HR, MR, LR, Safe)
%                                 Prob(1,m)= prob(beta,thresholds,C,m,DCPA,TCPA);
%                             end
%                             
%                             %Collision Risk Value
%                             CR(i,j)=  Prob(1,:) * risk_score(:,C);
%                         else
%                             CR(i,j)=0;
%                         end
                           
                    end
            
                end
            end
        end
    end
end
f_name = ['Result\','zone\','bivariate',t_indx,'.mat'];
save(f_name,'Matrx');
delete(['output_sync\','zone\','*.mat']);
delete(['output\','zone\','*.mat']);
end
end
     