clear
clc

a1{1}=[1,2];
a1{2}=2;

a2{1}=4;
a2{2}=5;
a2{3}=6;

nameList={a1,a2}

kk=1;
for ii = 1:length(nameList)
    for jj = 1:length(nameList{ii})
        a{kk}=nameList{ii}{jj};
        kk=kk+1;
    end
end

a
clearvars -except a