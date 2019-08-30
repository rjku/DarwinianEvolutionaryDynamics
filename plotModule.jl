
figure(1); clf()
plot(collect(1:NGEN),isingData.aveFitness); gcf()

figure(1.5); clf()
plot(collect(1:NGEN),isingData.growthFactor); gcf()

figure(1.6); clf()
plot(collect(1:NGEN),isingData.mutationNumber); gcf()

figure(2); clf()
matshow(aAves[1],cmap="Greys_r",vmin=-1,vmax=1); gcf()

figure(3); clf()
matshow(aAves[2],cmap="Greys_r",vmin=-1,vmax=1); gcf()

figure(4); clf()
# my_cmap=get_cmap("YlOrRd")
my_cmap=get_cmap("Greens")
my_cmap.set_under("Black")
my_cmap.set_over("White")

# matshow(JijMat,cmap=my_cmap,vmin=0,vmax=exp(maximum(isingPop.aGty[1].X))+.1);
matshow(JijMat,cmap=my_cmap,vmin=0,vmax=exp(maximum(myFancyX))+.1);
xticks(0:10:2SYSTEMSIZE,0:5:2SYSTEMSIZE÷2);
yticks(0:10:2SYSTEMSIZE,0:5:2SYSTEMSIZE÷2);
colorbar();
gcf()
