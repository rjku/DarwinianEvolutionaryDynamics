
using PyPlot, Statistics, LinearAlgebra

# *******************
# POPULATION DYNAMICS
# *******************

figure(1.1); clf(); title("Average Fitness"); plot( vcat([dataBtc.aveFitness for dataBtc in aIsingData]...) ); gcf()
figure(1.2); clf(); title("Growth Factor"); plot( vcat([dataBtc.growthFactor for dataBtc in aIsingData]...) ); gcf()
figure(1.3); clf(); title("Mutation Factor"); plot( vcat([dataBtc.mutationFactor for dataBtc in aIsingData]...) ); gcf()


# *******************
# GENOTYPE CHARACTERISTICS
# *******************

iBatch, iPop, iIsing = 1, 2, 1

GMat = Array{Float64}(undef, 2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L, 2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L)
@time mEvoFunc.showG!(aIsingData[iBatch].aLivingPop[iPop].aGty[1], GMat)

# figure(3); clf(); my_cmap=get_cmap("Greens"); my_cmap.set_under("Black"); my_cmap.set_over("White")
figure(3); clf(); my_cmap=get_cmap("YlOrRd"); my_cmap.set_under("Black"); my_cmap.set_over("White")
# figure(3); clf(); my_cmap=get_cmap("viridis"); my_cmap.set_under("Black"); my_cmap.set_over("White")
matshow(GMat,cmap=my_cmap,vmin=aIsingData[iBatch].aLivingPop[iPop].aGty[1].gbounds[1],vmax=aIsingData[iBatch].aLivingPop[iPop].aGty[1].gbounds[2]);
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1,
	1:Int32(aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].halfL):aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1,
	1:Int32(aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].halfL):aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
colorbar(); gcf()

# *******************
# SPIN STATISTICS
# *******************

const NSMPLSTAT, NMCSPSSTAT, NTRLSTAT = Int32(10^2), Int32(10^5), Int32(2)
isingDTMCprmStat = tDTMCprm( NSMPLSTAT, NMCSPSSTAT, NTRLSTAT )

aSpinAve = [ zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L
	) for i in isingEnv.IOidl ]
aSpinCov = [ zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2
	) for i in isingEnv.IOidl ]
sCor = Array{Float64}(undef, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2)

@time mEvoFunc.getSpinStat!(isingEnv, aIsingData[iBatch].aLivingPop[iPop].aGty[1], aSpinAve, aSpinCov, sCor, isingDTMCprmStat)

figure(2.1); clf(); matshow(aSpinAve[1],cmap="Greys_r",vmin=-1,vmax=1);
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L); gcf()

figure(2.2); clf(); matshow(aSpinAve[2],cmap="Greys_r",vmin=-1,vmax=1);
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L-1, 1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L); gcf()

figure(2.3); clf(); matshow(aSpinCov[1],cmap="viridis"); colorbar();
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2); gcf()

figure(2.4); clf(); matshow(aSpinCov[2],cmap="viridis"); colorbar();
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2); gcf()

figure(2.5); clf(); matshow(sCor,cmap="viridis",vmin=-1,vmax=1); colorbar();
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L:aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].L2); gcf()


# *******************
# Jij STATISTICS
# *******************

iBatch, iPop, iIsing = 4, 1, 1

JijAve = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dG)
JijCov = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dG, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dG)
JijCor = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dG, aIsingData[iBatch].aLivingPop[iPop].aGty[1].pMetaGty[1].dG)

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

iBatch, iPop, iIsing = 2, 5, 1

GAve = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].dG)
GCov = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].dG, aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].dG)
GCor = zeros(Float64, aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].dG, aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].dG)

mEvoFunc.getGStat!(aIsingData[iBatch].aLivingPop[iPop], GAve, GCov, GCor)

# figure(4.1); clf(); hist(vcat( [aIsingData[iBatch].aLivingPop[iPop].aGty[i].G for i in 1:aIsingData[iBatch].aLivingPop[iPop].pN[2]]... ),
# 	collect(-2:DELTAG:2),density=true,rwidth=.9); gcf()
# figure(4.1); clf(); hist(vcat( [aIsingData[iBatch].aLivingPop[iPop].aGty[i].G for i in 1:aIsingData[iBatch].aLivingPop[iPop].pN[2]]... ),
# 	aIsingData[iBatch].aLivingPop[iPop].aGty[1].g,density=true,rwidth=.9); gcf()
figure(4.1); clf(); hist(vcat( [aIsingData[iBatch].aLivingPop[iPop].aGty[i].G for i in 1:aIsingData[iBatch].aLivingPop[iPop].pN[2]]... ),
	collect(aIsingData[iBatch].aLivingPop[iPop].aGty[1].gbounds[1]:.2:aIsingData[iBatch].aLivingPop[iPop].aGty[1].gbounds[2]),
	density=true,rwidth=.9); gcf()
# figure(4.2); clf(); hist(vcat( [broadcast(x->10^(x),isingPop.aGty[i].G) for i in 1:isingPop.pN[2]]... ),4,density=true,rwidth=.9); gcf()

figure(5.1); clf(); matshow(GCor,cmap="viridis",vmin=-1,vmax=1); # title("Genotypic Variables Correlation");
xticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L2);
yticks(0:aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L2-1,
	1:aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L:2aIsingData[iBatch].aLivingPop[iPop].aGty[iIsing].pMetaGty[1].L2);
colorbar(); gcf()

evG = replace( e -> e < 10^-12 ? 0.0 : e, eigvals( Symmetric(GCov) ) )

figure(6.1); clf(); hist( evG, 40, density=true, log=true, rwidth=.9 ); gcf()
