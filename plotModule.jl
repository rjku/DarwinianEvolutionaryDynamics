
using PyPlot, Statistics, LinearAlgebra

iBatch, iPop, iIsing = 6, 1, 1

const NITSTAT, NSMPLSTAT, NTRLSTAT = Int32(10^3), Int32(10^3), Int32(2)
isingDTMCprmStat = tDTMCprm( NITSTAT, NSMPLSTAT, NTRLSTAT )

# *******************
# POPULATION DYNAMICS
# *******************

figure(1.1); clf(); title("Average Fitness"); plot(collect(1:NGEN),aIsingData[iBatch].aveFitness); gcf()
figure(1.2); clf(); title("Growth Factor"); plot(collect(1:NGEN),aIsingData[iBatch].growthFactor); gcf()
figure(1.3); clf(); title("Mutation Factor"); plot(collect(1:NGEN),aIsingData[iBatch].mutationFactor); gcf()


# *******************
# GENOTYPE CHARACTERISTICS
# *******************

JijMat = Array{Float64}(undef, 2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L, 2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L)
@time mEvoFunc.showJij!(aIsingData[iBatch].aLivingPop[iPop].aGty[1], JijMat)

# my_cmap=get_cmap("YlOrRd")
figure(3); clf(); my_cmap=get_cmap("Greens"); my_cmap.set_under("Black"); my_cmap.set_over("White")
matshow(JijMat,cmap=my_cmap,vmin=0,vmax=10^(maximum(aIsingData[iBatch].aLivingPop[iPop].aGty[1].X))+.1);
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1,
	1:Int32(aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].halfL):aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1,
	1:Int32(aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].halfL):aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
colorbar(); gcf()


# *******************
# SPIN STATISTICS
# *******************

aSpinAve = [ zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L) for i in isingEnv.idealInputOutput ]
aSpinCov = [ zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2) for i in isingEnv.idealInputOutput ]
sCor = Array{Float64}(undef, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2)

@time mEvoFunc.getSpinStat!(isingEnv, aIsingData[iBatch].aLivingPop[iPop].aGty[1], aSpinAve, aSpinCov, sCor, isingDTMCprmStat)

figure(2.1); clf(); matshow(aSpinAve[1],cmap="Greys_r",vmin=-1,vmax=1);
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L); gcf()

figure(2.2); clf(); matshow(aSpinAve[2],cmap="Greys_r",vmin=-1,vmax=1);
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L); gcf()

figure(2.3); clf(); matshow(aSpinCov[1],cmap="viridis"); colorbar(); gcf()
figure(2.4); clf(); matshow(aSpinCov[2],cmap="viridis"); colorbar(); gcf()
figure(2.5); clf(); matshow(sCor,cmap="viridis",vmin=-1,vmax=1); colorbar(); gcf()


# *******************
# Jij STATISTICS
# *******************

JijAve = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dX)
JijCov = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dX, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dX)
JijCor = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dX, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dX)

mEvoFunc.getJijStat!(aIsingData[iBatch].aLivingPop[iPop], JijAve, JijCov, JijCor)

figure(5.1); clf(); matshow(JijCor,cmap="viridis",vmin=-1,vmax=1);
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2);
colorbar(); gcf()

ev = eigvals(Symmetric(JijCov))

# figure(6.1); clf(); hist( map(l -> log(l), ev), 20, density=true, rwidth=.9); gcf()
figure(6.1); clf(); hist( ev, 20, density=true, rwidth=.9); gcf()

# *******************
# GENOTYPE STATISTICS
# *******************

XAve = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dX)
XCov = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dX, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dX)
XCor = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dX, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dX)

mEvoFunc.getXStat!(aIsingData[iBatch].aLivingPop[iPop], XAve, XCov, XCor)

figure(4.1); clf(); hist(vcat( [aIsingData[iBatch].aLivingPop[iPop].aGty[i].X for i in 1:aIsingData[iBatch].aLivingPop[iPop].pN[2]]... ),
	collect(-4:DELTAX:4),density=true,rwidth=.9); gcf()
# figure(4.2); clf(); hist(vcat( [broadcast(x->10^(x),isingPop.aGty[i].X) for i in 1:isingPop.pN[2]]... ),4,density=true,rwidth=.9); gcf()

figure(5.1); clf(); matshow(XCor,cmap="viridis",vmin=-1,vmax=1);
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L2);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L2);
colorbar(); gcf()

evX = eigvals(Symmetric(XCov))

figure(6.1); clf(); hist( evX, 30, density=true, rwidth=.9 ); gcf()
