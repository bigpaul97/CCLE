function ave = aveHiC(sampleIndex,fileName)
% assume sampleIndex has at least 2 elements

load(append(fileName,'_',num2str(sampleIndex(1)),'.mat'),'hmap')

hmapT = hmap;

for i = sampleIndex(2:end)
    load(append(fileName,'_',num2str(i),'.mat'),'hmap')
    hmapT = hmapT + hmap;
end

ave = hmapT./length(sampleIndex);