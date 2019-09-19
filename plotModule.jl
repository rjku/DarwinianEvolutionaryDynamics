
using PyPlot

figure(1.1); clf()
title("Average Fitness")
plot(collect(1:NGEN),aIsingData[end].aveFitness); gcf()

figure(1.2); clf()
title("Growth Factor")
plot(collect(1:NGEN),aIsingData[end].growthFactor); gcf()

figure(1.3); clf()
title("Mutation Factor")
plot(collect(1:NGEN),aIsingData[end].mutationFactor); gcf()

figure(2.1); clf()
matshow(aAves[1],cmap="Greys_r",vmin=-1,vmax=1); gcf()

figure(2.2); clf()
matshow(aAves[2],cmap="Greys_r",vmin=-1,vmax=1); gcf()

# using PyPlot
figure(3); clf()
# my_cmap=get_cmap("YlOrRd")
my_cmap=get_cmap("Greens")
my_cmap.set_under("Black")
my_cmap.set_over("White")

matshow(JijMat,cmap=my_cmap,vmin=0,vmax=exp(maximum(isingPop.aGty[1].X))+.1);
# matshow(JijMat,cmap=my_cmap,vmin=0,vmax=exp(maximum(isingPopClone.aGty[1].X))+.1);
# matshow(JijMat,cmap=my_cmap,vmin=0,vmax=exp(maximum(myFancyX))+.1);
xticks(0:10:2isingPop.aGty[1].pMetaGty[1].L,0:5:isingPop.aGty[1].pMetaGty[1].L);
yticks(0:10:2isingPop.aGty[1].pMetaGty[1].L,0:5:isingPop.aGty[1].pMetaGty[1].L);
colorbar();
gcf()

figure(4.1); clf()
hist(isingPop.aGty[1].X,20,rwidth=.85)
gcf()

figure(4.2); clf()
hist(broadcast(x->exp(x),isingPop.aGty[1].X),8,rwidth=.85)
gcf()
