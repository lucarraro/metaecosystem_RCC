rm(list=ls())

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


library(deSolve)
library(OCNet)
library(R.matlab)


cat("Loading OCN...\n")
if (!file.exists("support/OCN_750.rda")){
  set.seed(1)
  OCN <- create_OCN(750,750, nIter = 80 * 1000 * 1000, # WARNING! This takes weeks to run
                    typeInitialState = "V",    
                    initialNoCoolingPhase = 0.05,    
                    coolingRate = 0.25,
                    displayUpdates = 2,
                    nUpdates = 1000,
                    saveEnergy = TRUE)
  OCN$energy <- NULL
  cellsize <- 100
  OCN$cellsize <- cellsize
  OCN$FD$A <- OCN$FD$A*cellsize^2
  OCN$FD$X <- cellsize/2 + (OCN$FD$X - min(OCN$FD$X))*cellsize
  OCN$FD$Y <- cellsize/2 + (OCN$FD$Y - min(OCN$FD$Y))*cellsize
  
  OCN <- landscape_OCN(OCN,slope0 = 0.0035,displayUpdates = 2)
  
  OCN <- aggregate_OCN(OCN, thrA=1e6, maxReachLength=2000) 
  
  OCN <- rivergeometry_OCN(OCN, widthMax = 24, depthMax = 5, velocityMax = 1.25)
  
  OCN$AG$leng[OCN$AG$outlet] <- 500
  
  save(file="support/OCN_750.rda",OCN)
  
} else {
  load("support/OCN_750.rda")
}

Q <- OCN$AG$width * OCN$AG$depth * OCN$AG$velocity
V <- OCN$AG$width * OCN$AG$depth * OCN$AG$leng

writeMat("support/OCNparam.mat",downNode=OCN$AG$downNode,AS=OCN$SC$A,B=OCN$AG$width,L=OCN$AG$leng,
         Q=Q,V=V,A=OCN$AG$A,streamOrder=OCN$AG$streamOrder,Z=OCN$FD$Z,subcatchMat=OCN$FD$toSC, depth=OCN$AG$depth)


# Fig 2a
cat("Producing Fig. 2a...\n")
pdf('Fig2a.pdf',width=18,height=12,paper="special")
draw_simple_OCN(OCN, thrADraw = OCN$thrA)

nodeAG_vec <- c(958,1601,1913)
for (nodeAG in nodeAG_vec){
up_nodes <- OCN$AG$upstream[[nodeAG]]
nodes <- numeric(0)
for (i in 1:length(up_nodes)){
  nodes <- c(nodes,OCN$SC$toFD[[up_nodes[i]]])
}

subset_X <- OCN$FD$X[nodes]
subset_Y <- OCN$FD$Y[nodes]
X_mesh <- seq(min(subset_X)-OCN$cellsize,max(subset_X)+OCN$cellsize,OCN$cellsize)
Y_mesh <- seq(min(subset_Y)-OCN$cellsize,max(subset_Y)+OCN$cellsize,OCN$cellsize)
mesh <- matrix(data=0,nrow=length(Y_mesh),ncol=length(X_mesh))
for (i in 1:length(subset_X)){
  ind_x <- which(X_mesh==subset_X[i])
  ind_y <- which(Y_mesh==subset_Y[i])
  mesh[ind_y,ind_x] <- 1
}
count <- contourLines(X_mesh,Y_mesh,t(mesh),levels=1)
X_contour <- vector("list", length = length(count))
Y_contour <- vector("list", length = length(count))
for (k in 1:length(count)){
  X_contour[[k]] <- count[[k]]$x
  Y_contour[[k]] <- count[[k]]$y
} 

  for (k in 1:length(X_contour)){
    lines(X_contour[[k]],Y_contour[[k]],lwd=2,col="black")
  }
}
dev.off()

# Fig. 3 a--d
cat("Producing Fig. 3 (maps)...\n")

if (file.exists("support/resultsR.mat")){
  res <- readMat("support/resultsR.mat")
} else {
  stop("resultsR.mat missing. You need to run MAIN_mixed.m")
}


indG <- seq(2,10*OCN$AG$nNodes,by=10); indS <- seq(3,10*OCN$AG$nNodes,by=10); 
indC <- seq(4,10*OCN$AG$nNodes,by=10); indF <- seq(5,10*OCN$AG$nNodes,by=10);

pdf('Fig3a--d.pdf',width=18,height=12,paper="special")
par(mfrow=c(2,2))
cat("             Grazers\n")
draw_thematic_OCN(res$X[indG],OCN,colLevels=c(0,2e-6,1000)); title("Grazers")
cat("             Shredders\n")
draw_thematic_OCN(res$X[indS],OCN,colLevels=c(0,2e-6,1000)); title("Shredders")
cat("             Collectors\n")
draw_thematic_OCN(res$X[indC],OCN,colLevels=c(0,2e-6,1000)); title("Collectors")
cat("             Filter feeders\n")
draw_thematic_OCN(res$X[indF],OCN,colLevels=c(0,2e-6,1000)); title("Filter feeders")
dev.off()
