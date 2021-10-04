function makegrid(gridindex,dim)
hold on
for i = 1:length(gridindex)
    v = gridindex(i);
    base = [0 dim];
    plot([v v],base,'color', [.5 .5 .5]);
    plot(base,[v v],'color', [.5 .5 .5]);
end
hold off