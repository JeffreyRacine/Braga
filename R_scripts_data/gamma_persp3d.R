load("~/tex/Rennes_2022/work/sandbox_persp/skip_M_400_sigma_0.50/mc_online.RData")

library(rgl)
options(rgl.printRglwidget = TRUE)

open3d()
mfrow3d(2,2, sharedMouse = TRUE)
par3d(windowRect=c(0,0,800,800),cex=0.75)

num.colors <- 1000

z <- fda.la.out$Gamma.hat
colorlut <- topo.colors(num.colors)
col <- colorlut[ (num.colors-1)*(z-min(z))/(max(z)-min(z)) + 1 ]

x1.seq <- x2.seq <- grid.t0
persp3d(x=x1.seq,y=x2.seq,z=z,
        xlab="s",ylab="t",zlab="Gamma(s,t)",
        ticktype="detailed",
        border="red",
        color=col,
        alpha=.7,
        main="Discard, No Correction",
        sub=paste("RMSE = ",format( RMSE.Gamma[m,mc],digits=4),
                  sep=""))
grid3d(c("x", "y+", "z"))


z <- fda.la.out$Gamma.hat.flattop.correction
colorlut <- topo.colors(num.colors)
col <- colorlut[ (num.colors-1)*(z-min(z))/(max(z)-min(z)) + 1 ]

x1.seq <- x2.seq <- grid.t0
persp3d(x=x1.seq,y=x2.seq,z=z,
        xlab="s",ylab="t",zlab="Gamma(s,t)",
        ticktype="detailed",
        border="red",
        color=col,
        alpha=.7,
        main="Discard",
        sub=paste("RMSE = ",format( RMSE.Gamma.flattop[m,mc],digits=4),
                  sep=""))
grid3d(c("x", "y+", "z"))

load("~/tex/Rennes_2022/work/sandbox_persp/lp_M_400_sigma_0.50/mc_online.RData")

z <- fda.la.out$Gamma.hat.flattop.correction
colorlut <- topo.colors(num.colors)
col <- colorlut[ (num.colors-1)*(z-min(z))/(max(z)-min(z)) + 1 ]

x1.seq <- x2.seq <- grid.t0
persp3d(x=x1.seq,y=x2.seq,z=z,
        xlab="s",ylab="t",zlab="Gamma(s,t)",
        ticktype="detailed",
        border="red",
        color=col,
        alpha=.7,
        main="Reconstruct",
        sub=paste("RMSE = ",format( RMSE.Gamma.flattop[m,mc],digits=4),
                  sep=""))
grid3d(c("x", "y+", "z"))

z <- Gamma.true
colorlut <- topo.colors(num.colors)
col <- colorlut[ (num.colors-1)*(z-min(z))/(max(z)-min(z)) + 1 ]

persp3d(x=x1.seq,y=x2.seq,z=z,
        xlab="s",ylab="t",zlab="Gamma(s,t)",
        ticktype="detailed",
        border="red",
        color=col,
        alpha=.7,
        main="True")
grid3d(c("x", "y+", "z"))
