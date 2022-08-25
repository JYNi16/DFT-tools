%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dig_curve_data version3.0
% Author: Xiao, Ruichu
% Address: Institute of Solid State Physics, Chinese Academy of Sciences
% Email: xiaoruichun@foxmail.com
% Blog:  http://blog.sciencenet.cn/u/XRC0808087
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
format long
 
%%%%%%%%%%%%%%%%设定需要处理的能带序号
nband_start=2;
nband_end=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   读取BXSF文档，取得kprim efermi  nx/ny/nz
[filename,pathname]=uigetfile('*BXSF');
infile=strcat(pathname,filename);
fid=fopen(infile,'r');
while  ~feof(fid) 
    str = fgetl(fid);   %读取一行   
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%读取BXSF文档，获得EIG数据写入临时文件
fid=fopen(infile,'r');
fid1=fopen('data_tmp.txt','a');
while  ~feof(fid)                   %前18行为文件有用信息
     str = fgetl(fid);   %读取一行
      if ~isempty(strfind(str,'BAND:')) 
          break;
      end
end
while  ~feof(fid) 
   str = fgetl(fid);   %读取一行       
      if  isempty(strfind(str,'BAND')) &   isempty(strfind(str,'END'))  %%(str(1)~='B' & str(2)~='E') 
        fprintf(fid1,'%s \n',str); %将一行数据写入文件  
      end
end
fclose(fid);
fclose(fid1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%将能量数据存储到矩阵 EIG中
fid1=fopen('data_tmp.txt','r');
data_bxsf=fscanf(fid1,'%f');
i=1;
for m=1:nband
    for i_x = 1:nx
        for i_y = 1:ny
            for i_z =1:nz    
                EIG(i_x,i_y,i_z,m)=data_bxsf(i)-E_fermi;    %已将费米能级设置为零参考,并转为eV
                i=i+1;               
            end 
        end
    end
end
fclose(fid1);
delete('data_tmp.txt')   
%%%%%尚未编写完成
long=nx*ny*nz;
for i=1:nband
   range(i,:)=[ min(min(min(EIG(:,:,:,i))))  ,max(max(max(EIG(:,:,:,i))))];
   %data_bxsf(1:long);
end
%%%%%%%%%
n_inter=input('输入内插倍数(建议2-5倍):  ');
NX=(nx-1)*n_inter+1;
NY=(ny-1)*n_inter+1;
NZ=(nz-1)*n_inter+1;

x=linspace(0,1,nx);
y=linspace(0,1,ny);
z=linspace(0,1,nz);
X=linspace(0,1,NX);
Y=linspace(0,1,NY);
Z=linspace(0,1,NZ);
[xx,yy,zz]=meshgrid(y,x,z);
[XX,YY,ZZ]=meshgrid(Y,X,Z);

%E_inter(:,:,:,m)=interp3(EIG(:,:,:,m),1,'cubic');
for i=1:nband
E_inter(:,:,:,i)=interp3(xx,yy,zz,EIG(:,:,:,i),XX,YY,ZZ,'spline');
end
%E_inter=imresize(EIG,2,'nearest');
%[NX,NY,NZ,m]=size(E_inter);
for m=nband_start:nband_end
    j=1; data_xsf=[];
    for i_z =1:NZ
      for i_y = 1:NY
        for i_x = 1:NX
                data_xsf(j)=E_inter(i_x,i_y,i_z,m);    %转换为xsf文件DATAGRIES输出数据的顺序：
                j=j+1;               
            end 
        end
    end


outfile=strcat(pathname,filename,'_band_',char(num2str(m)),'.xsf');
fid2=fopen(outfile,'w');
%%%%写入xsf文件数据
fprintf(fid2,'BEGIN_BLOCK_DATAGRID_3D \n');
fprintf(fid2,'bxsf2xsf_created_by_RC_Xiao \n');
fprintf(fid2,'BEGIN_DATAGRID_3D \n');
fprintf(fid2,'%d %d  %d  \n',NX,NY,NZ);
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

