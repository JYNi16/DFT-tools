%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dig_curve_data version3.0
% Author: Xiao, Ruichu
% Address: Institute of Solid State Physics, Chinese Academy of Sciences
% Email: xiaoruichun@foxmail.com
% Blog:  http://blog.sciencenet.cn/u/XRC0808087
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
format long
%%%%%%%%%%%%%%%%�趨��Ҫ������ܴ����
nband_start=1;
nband_end=13;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ��ȡBXSF�ĵ���ȡ��kprim efermi  nx/ny/nz
[filename,pathname]=uigetfile('*BXSF');
infile=strcat(pathname,filename);
fid=fopen(infile,'r');
while  ~feof(fid) 
    str = fgetl(fid);   %��ȡһ��   
    S = regexp(str, '\s+', 'split');
    if ~isempty(str)
        if ~isempty(strfind(str,'Fermi'))     
            E_fermi=str2num(char(S(4)));    
        end    
        
    if ~isempty(strfind(str,'BANDGRID_3D_BANDS'))   
         str = fgetl(fid);  
         S = regexp(str, '\s+', 'split');
         nband=str2num(char(S(2))); 
         str = fgetl(fid);     
         S = regexp(str, '\s+', 'split');
         nx=str2num(char(S(2))); 
         ny=str2num(char(S(3))); 
         nz=str2num(char(S(4))); 
         str = fgetl(fid);  
         str = fgetl(fid);   
         S = regexp(str, '\s+', 'split');
         kx=[str2num(char(S(2))), str2num(char(S(3))), str2num(char(S(4)))];
         str = fgetl(fid);  
         S = regexp(str, '\s+', 'split');
         ky=[str2num(char(S(2))), str2num(char(S(3))), str2num(char(S(4)))];
         str = fgetl(fid);  
         S = regexp(str, '\s+', 'split');
         kz=[str2num(char(S(2))), str2num(char(S(3))), str2num(char(S(4)))];
         break;
    end    
    end
end

fclose(fid);
ngrid=[nx ny nz];
fid1=fopen('data_tmp.txt','w');
fclose(fid1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%��ȡBXSF�ĵ������EIG����д����ʱ�ļ�
fid=fopen(infile,'r');
fid1=fopen('data_tmp.txt','a');
while  ~feof(fid)                   %ǰ18��Ϊ�ļ�������Ϣ
     str = fgetl(fid);   %��ȡһ��
      if ~isempty(strfind(str,'BAND:')) 
          break;
      end
end
while  ~feof(fid) 
   str = fgetl(fid);   %��ȡһ��       
      if  isempty(strfind(str,'BAND')) &   isempty(strfind(str,'END'))  %%(str(1)~='B' & str(2)~='E') 
        fprintf(fid1,'%s \n',str); %��һ������д���ļ�  
      end
end
fclose(fid);
fclose(fid1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���������ݴ洢������ EIG��
fid1=fopen('data_tmp.txt','r');
data_bxsf=fscanf(fid1,'%f');
i=1;
for m=1:nband
    for i_x = 1:nx
        for i_y = 1:ny
            for i_z =1:nz    
                EIG(i_x,i_y,i_z,m)=data_bxsf(i)-E_fermi;    %�ѽ������ܼ�����Ϊ��ο�,��תΪeV
                i=i+1;               
            end 
        end
    end
end
fclose(fid1);
delete('data_tmp.txt')   
%%%%%��δ��д��ɣ�����Ӱ��ʹ��
long=nx*ny*nz;
for i=1:nband
   range(i,:)=[ min(min(min(EIG(:,:,:,i))))  ,max(max(max(EIG(:,:,:,i))))];
   %data_bxsf(1:long);
end
%%%%%%%%%

for m=nband_start:nband_end
    j=1; data_xsf=[];
    for i_z =1:nz
      for i_y = 1:ny
        for i_x = 1:nx
                data_xsf(j)=EIG(i_x,i_y,i_z,m);    %ת��Ϊxsf�ļ�DATAGRIES������ݵ�˳��
                j=j+1;               
            end 
        end
    end


outfile=strcat(pathname,filename,'_band_',char(num2str(m)),'.xsf');
fid2=fopen(outfile,'w');
%%%%д��xsf�ļ�����
fprintf(fid2,'BEGIN_BLOCK_DATAGRID_3D \n');
fprintf(fid2,'bxsf2xsf_created_by_RC_Xiao \n');
fprintf(fid2,'BEGIN_DATAGRID_3D \n');
fprintf(fid2,'%d %d  %d  \n',ngrid);
fprintf(fid2,'0.000000  0.000000  0.000000 \n');
fprintf(fid2,'%f  %f  %f \n',kx);
fprintf(fid2,'%f  %f  %f \n',ky);
fprintf(fid2,'%f  %f  %f \n',kz);
fprintf(fid2,'%f  %f   %f   %f   %f   %f \n',data_xsf);
fprintf(fid2,'\n');
fprintf(fid2,'END_DATAGRID_3D \n');
fprintf(fid2,'END_BLOCK_DATAGRID_3D \n');
end

status='finished!'

