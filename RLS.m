yitam=0.1;
yitad=0.1;
%function [result,result2]=mymain(yitad,yitam)
%predict disease-related miRNA based on RLS
%A: Binary relations between disease and miRNA, 1st column:disease, 2nd column:miRNA
A=textread('knowndiseasemirnainteraction.txt');
% nd:the number of diseases
% nm:the number of miRNAs
% pp:the number of known diseae-miRNA associations
nd=max(A(:,2));
nm=max(A(:,1));
[pp,qq]=size(A);
%interaction: adajency matrix for the disease-miRNA association network
%interaction(i,j)=1 means miRNA j is related to disease i
rank=ones(5430,1);
rankl=ones(5430,1);
interaction=zeros(495,383);
for i=1:pp
    interaction(A(i,1),A(i,2))=1;
end
interaction=interaction';
id=0;

               %clc
           %fprintf('ii=%d,jj=%d\n',ii,jj);
               id=id+1;
           fprintf('id=%d\n',id);


%load diseasesimilarity;
%kd=diseasesimilarity;
kd=textread('K.txt');
%functionalsimilarity=textread('adajency matrix.txt');
%km=functionalsimilarity;
km=textread('K2.txt');
%obatian predictive result in the disease space 

F1=km*pinv(km+yitam*eye(nm))*interaction';
F1=F1';
%obatian predictive result in the miRNA space 
F2=kd*pinv(kd+yitad*eye(nd))*interaction;
%obatian predictive result in the combined space 
F=F1*0.9+F2*0.1;

for i1=1:nd
    for j=1:nm
        if interaction(i1,j)==1
           F(i1,j)=-100000;
        end
    end
end

finalprediction=[];
nm=495;
nd=383;
for i=1:nm
    finalprediction=[finalprediction;F(:,i)];
end

finalpredictiondisease=[];
for i=1:nm
    finalpredictiondisease=[finalpredictiondisease,1:nd];
end
finalprediction(:,2)=finalpredictiondisease';
finalpredictionmirna=[];
for i=1:nm
    finalpredictionmirna=[finalpredictionmirna;i.*ones(nd,1)];
end
finalprediction(:,3)=finalpredictionmirna;
%sortrows mean ranking according to given column,here 1st column, -1 means descend rank
result=sortrows(finalprediction,-1);
for i=1:383
id=1;
for j=1:383*495
if result(j,2)==i
final1(i,id)=result(j,1);
id=id+1;
end
end
end
for i=1:383
id=1;
for j=1:383*495
if result(j,2)==i
final2(i,id)=result(j,3);
id=id+1;
end
end
end
id=1;
for i=1:495*383
total(i,3)=num2cell(result(i,1));
total(i,1)=num2cell(result(i,2));
total(i,2)=num2cell(result(i,3));
%total(i,1)=diseasename(result(i,2));
%total(i,2)=miRNAname(result(i,3));
end
