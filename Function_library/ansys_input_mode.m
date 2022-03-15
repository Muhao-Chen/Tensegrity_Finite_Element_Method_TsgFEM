function[]= ansys_input_mode(N,C,A,t,b,index_s,name)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%this function out put the code used in ANSYS from matlab
%% predefinition for output
[ne,nn]=size(C);
N_xyz_nn=N';
for i=1:ne
    El_123_nn(i,:)=find(C(i,:));
end
prestress=t./A;
%% output for ansys
fid11=fopen(name,'w');
%%
fprintf(fid11,'!输出节点坐标，拓扑关系，截面积，预应力\n!材料特性，弹塑性需要自己更改\n\n');
fprintf(fid11,'finish\n/clear \n/filename,tower  \n/title,the analysis of tower  \n!单位m，N，Pa，s\n\n');

%%  material properties
fprintf(fid11,'/prep7\n');
fprintf(fid11,'!定义单元类型 \n et,1,link180 \n \n');
fprintf(fid11,'!定义基本参数(数据) \n es=76000e6 \n eg=2.06e11 \nfd=2.15e08 \nfy=2.35e08 \n ft=1.57e09 \n ratio1=0.85 \nratio2=0.70\n');
fprintf(fid11,'!定义材料特性(数据) \n mp,ex,1,es	!定义弹性模量 \n mp,prxy,1,0.3	!定义主泊松比\n mp,dens,1,7870	!定义质量密度\nmp,alpx,1,6.5e-6	!定义线膨胀系数\n\n');
fprintf(fid11,'!定义材料特性(数据) \n mp,ex,2,eg	!定义弹性模量 \n mp,prxy,2,0.3	!定义主泊松比\n mp,dens,2,7870	!定义质量密度\nmp,alpx,2,6.5e-6	!定义线膨胀系数\n\n');

%% nodal coordinates,connectivity

for i=1:nn
    fprintf(fid11,'K,%d,%17.15f,%17.15f,%17.15f  !节点坐标\n',i,N_xyz_nn(i,:));
end
fprintf(fid11,'\n');
for i=1:ne
    fprintf(fid11,'L,%4d,%4d  !创建线\n',El_123_nn(i,:));
end
fprintf(fid11,'\n');
%% area
fprintf(fid11,'*dim,area,,%d\n',ne);
for i=1:ne
    fprintf(fid11,'area(%d)=%4d !截面积\n',i,A(i));
end
fprintf(fid11,'\n');
%  fprintf(fid11,'*DO,J,1,%d\n sectype,J,link  !定义截面类型为杆\nsecdata,area(J)	!定义杆截面几何数据\nseccontrol,,%d	!定义只拉不压（1）可拉可压（0）\n*ENDDO\n',ne,max(i==index_s));

for i=1:ne
    fprintf(fid11,'sectype,%d,link  !定义截面类型为杆\nsecdata,area(%d)   !定义杆截面几何数据\nseccontrol,,%d       !定义只拉不压（1）可拉可压（0）\n',i,i,max(i==index_s));
end
fprintf(fid11,'\n');
%% 定义杆件属性
for i=1:ne
    fprintf(fid11,'!定义单元属性\nlsel,s,,,%d  !选择杆件\nlatt,%d,,1,,,,%d  !定义杆件截面积\n',i,2-max(i==index_s),i);
end
fprintf(fid11,'\n');
%% prestress
fprintf(fid11,'*dim,prestress,,%d\n',ne);
for i=1:ne
    fprintf(fid11,' prestress(%d)=%4f  !预应力\n',i,prestress(i));
end
fprintf(fid11,'\n');
%% mesh
fprintf(fid11,'!划分单元 \n LSEL,ALL \n LESIZE,ALL,,,1\nLMESH,ALL\nfinish\n');
fprintf(fid11,'\n');
%% solve
fprintf(fid11,'!第一次求解（自平衡）\n/SOLU\nANTYPE,0 \nNLGEO!考虑大变形 \nSSTIF,ON	!考虑应力刚化 \nNSUBST,100	!设置荷载步的子步数 \nAUTOTS,ON	!自动荷载步 \n  OUTRES,ALL,ALL 	!输出控制 \n');
fprintf(fid11,'\n');
%% boundary constraints

for i=1:numel(b)
    switch b(i)-3*(ceil(b(i)/3)-1)
        case 1
            fprintf(fid11,'DK,%d,UX\n',ceil(b(i)/3));
        case 2
            fprintf(fid11,'DK,%d,UY\n',ceil(b(i)/3));
        case 3
            fprintf(fid11,'DK,%d,UZ\n',ceil(b(i)/3));
    end
end
fprintf(fid11,'\n');


%% prestress and solve
fprintf(fid11,'*DO,J,1,%d	!数据填入预应力数组\n	INISTATE,DEFINE,J,,,,PRESTRESS(J)\n*ENDDO\n',ne);
fprintf(fid11,'\n');
fprintf(fid11,'ALLSEL,ALL\nSOLVE\nFINISH\n');

fprintf(fid11,'\n');
%% post1
fprintf(fid11,'!后处理\n/POST1\nPLDISP !显示变形后的图形\nALLSEL,ALL  !显示内力云图\nETABLE,MFORCE,SMISC,1\nPLLS,MFORCE,MFORCE,0.4\n',ne);

%% 模态分析
%  fprintf(fid11,'!模态分析\n/SOLU\nANTYPE,2\nPRTRES,ON		!打开预应力效应\nMODOPT,LANB,18\nPLLS,MFORCE,MFORCE,0.4\n',ne);
%
% MXPAND,18,,,1
% /STATUS,SOLU	!查看求解SOLUTION OPTION
% SOLVE
%
% !在求解层可获得各阶模态频率、参与系数、模态系数、阻尼比等参数
% *DIM,FI,,18
% *DIM,PFI,,18
% *DIM,MCI,,18
% *DIM,DAI,,18
% *DO,I,1,18
% *GET,FI(I),MODE,I,FREQ
% *GET,PFI(I),MODE,I,PFACT
% *GET,MCI(I),MODE,I,MCOEF
% *GET,DAI(I),MODE,I,DAMP
% *ENDDO
% FINISH
%
% /POST1
% !SET,1,2
% !PLDISP,1
% !PLNSOL,U,Y
% !ETABLE,MFORCE,SMISC,1
% !PLLS,MFORCE,MFORCE,0.4
%
% SET,LIST
%
% !*CFOPEN,LOAD_DATA,TXT
%
% *MWRITE,FI,FI,TXT
% (F20.10)
%
%
% !*CFCLOS,LOAD_DATA
fclose(fid11);
end

