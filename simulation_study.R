###############
# SIMU 100 REP
###############

library(mvtnorm)
library(doParallel)
library(BNPmix)
library(mcclust.ext)
library(reshape)
library(viridis)
library(ggplot2)
library(ggpubr)

# Get upper triangle of the matrix
get_upper_tri <- function(mat) {
  mat[lower.tri(mat)] <- NA
  return(mat)
}

###############

set.seed(42)
num  <- c(100, 300, 1000)
nrep <- 1
cost <- c(1/5, 1/2, 1, 2, 5)
grid <- expand.grid(seq(-6.5,6.5,length.out = 40),seq(-6.5,6.5,length.out = 40))
ngrid <- nrow(grid)

nrep <- 100
using_list <- list()
k <- 1

for(l in 1:nrep){
  for(i in 1:length(num)){
    
    temp <- rbind(rmvnorm(num[i]/2, c(-2,-2), matrix(c(1, .85, .85, 1), ncol = 2)), rmvnorm(num[i]/2, c(2,2), diag(1,2)))
    
    for(j in 1:length(cost)){
      
      using_list[[k]] <- list()
      using_list[[k]]$grid <- as.matrix(grid) %*% diag(cost[j], 2)
      using_list[[k]]$data <- as.matrix(temp) %*% diag(cost[j], 2)
      using_list[[k]]$cost <- diag(cost[j], 2)
      using_list[[k]]$num <- num[i]
      using_list[[k]]$dataorig <- temp
      using_list[[k]]$mean <- colMeans(temp)
      using_list[[k]]$var <- var(temp)
      k <- k+1 
      
    }
  }
}
dens <- apply(grid, 1, function(x) dmvt(x, df = 2, log = F))

###################################

no_cores <- 6
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores)

a <- Sys.time()
result <- foreach(i = 1:1500, .packages = c("BNPmix", "mcclust.ext")) %dopar% {
  mod <- DPmixMulti(nsim = 7500,
                    nburn = 2500,
                    napprox = 100,
                    data = as.matrix(using_list[[i]]$data),
                    grid = as.matrix(using_list[[i]]$grid),
                    mu_start = c(0,0),
                    Lambda_start = diag(1,2),
                    theta = 1, 
                    m0 = c(0,0),
                    B0 = diag(1,2),
                    nu0 = 4,
                    sigma = diag(1,2),
                    b1 = 4,
                    B1 = diag(15,2),
                    m1 = c(0,0),
                    k1 = 1,
                    nupd = 200, 
                    fix = TRUE, 
                    dep = T)
  
  NGROUP <- mean(apply(mod[[2]], 1, function(x) length(table(x))))
  psm      <- comp.psm(mod[[2]] + 1)
  NGROUPvi <- length(table(minVI(psm, mod[[2]], method=("avg"), include.greedy=FALSE)$cl))
  
  list(mod[[1]], NGROUP, NGROUPvi)
}
stopCluster(cl) 
b <- Sys.time()
b - a

##### post simulation analysis 

ten100  <- array(0, dim = c(5, 5, 100))
ten300  <- array(0, dim = c(5, 5, 100))
ten1000 <- array(0, dim = c(5, 5, 100))
NGR   <- array(0, dim = c(5, 3, 100))
NGRVI <- array(0, dim = c(5, 3, 100))
grid <- as.matrix(grid)
for(i in 1:100){
  for(j in 1:5){
    for(k in 1:5){
      ten100[j,k,i] <- sum(abs(result[[(((i - 1) * 15 ) + j)]][[1]] * det(using_list[[(((i - 1) * 15 ) + j)]]$cost) - 
                                 result[[(((i - 1) * 15 ) + k)]][[1]] * det(using_list[[(((i - 1) * 15 ) + k)]]$cost))) * 
                                    ((unique(grid[,1])[3] - unique(grid[,1])[2]) * (unique(grid[,2])[3] - unique(grid[,2])[2])) 
      ten300[j,k,i] <- sum(abs(result[[(((i - 1) * 15 ) + j + 5)]][[1]] * det(using_list[[(((i - 1) * 15 ) + j + 5)]]$cost) - 
                                 result[[(((i - 1) * 15 ) + k + 5)]][[1]] * det(using_list[[(((i - 1) * 15 ) + k + 5)]]$cost))) * 
                                    ((unique(grid[,1])[3] - unique(grid[,1])[2]) * (unique(grid[,2])[3] - unique(grid[,2])[2])) 
      ten1000[j,k,i] <- sum(abs(result[[(((i - 1) * 15 ) + j + 10)]][[1]] * det(using_list[[(((i - 1) * 15 ) + j + 10)]]$cost) - 
                                 result[[(((i - 1) * 15 ) + k + 10)]][[1]] * det(using_list[[(((i - 1) * 15 ) + k + 10)]]$cost))) * 
                                    ((unique(grid[,1])[3] - unique(grid[,1])[2]) * (unique(grid[,2])[3] - unique(grid[,2])[2])) 
    }
    NGR[j,1,i] <- result[[(((i - 1) * 15 ) + j)]][[2]] 
    NGR[j,2,i] <- result[[(((i - 1) * 15 ) + j + 5)]][[2]]
    NGR[j,3,i] <- result[[(((i - 1) * 15 ) + j + 10)]][[2]]
    
    NGRVI[j,1,i] <- result[[(((i - 1) * 15 ) + j)]][[3]] 
    NGRVI[j,2,i] <- result[[(((i - 1) * 15 ) + j + 5)]][[3]]
    NGRVI[j,3,i] <- result[[(((i - 1) * 15 ) + j + 10)]][[3]]
  }
}
round(t(apply(NGR, c(1,2), mean)), digits = 2)
t(apply(NGRVI, c(1,2), mean))

# Plots

upper_tri_ten100 <- get_upper_tri(apply(ten100, c(1,2), mean))
melted_ten100 <- melt(upper_tri_ten100, na.rm = TRUE)
ggheatmap_ten100 <- ggplot(melted_ten100, aes(X2, X1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "D", name = "", direction=-1, na.value = "white", limits = c(0,0.81), guide = F) +
  coord_fixed() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size=24),
    plot.subtitle = element_text(hjust = 0.5, size=11)
  ) +
  guides(fill = guide_colorbar(
    barwidth = 1,
    barheight = 7,
    title.hjust = 0.5
  )) +
  scale_x_discrete(limit = c("1/5", "1/2", "1", "2", "5"))+ 
  scale_y_discrete(limit = c("1/5", "1/2", "1", "2", "5"))
ggheatmap_ten100

upper_tri_ten300 <- get_upper_tri(apply(ten300, c(1,2), mean))
melted_ten300 <- melt(upper_tri_ten300, na.rm = TRUE)
ggheatmap_ten300 <- ggplot(melted_ten300, aes(X2, X1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "D", name = "", direction=-1, na.value = "white", limits = c(0,0.81), guide = FALSE) +
  coord_fixed() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size=24),
    plot.subtitle = element_text(hjust = 0.5, size=11)
  ) +
  guides(fill = guide_colorbar(
    barwidth = 1,
    barheight = 7,
    title.position = "top",
    title.hjust = 0.5
  )) +
  scale_x_discrete(limit = c("1/5", "1/2", "1", "2", "5"))+ 
  scale_y_discrete(limit = c("1/5", "1/2", "1", "2", "5"))
ggheatmap_ten300

upper_tri_ten1000 <- get_upper_tri(apply(ten1000, c(1,2), mean))
melted_ten1000 <- melt(upper_tri_ten1000, na.rm = TRUE)
ggheatmap_ten1000 <- ggplot(melted_ten1000, aes(X2, X1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "D", name = "", direction=-1, na.value = "white", limits = c(0,0.81)) +
  coord_fixed() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size=24),
    plot.subtitle = element_text(hjust = 0.5, size=11)
  ) +
  scale_x_discrete(limit = c("1/5", "1/2", "1", "2", "5"))+ 
  scale_y_discrete(limit = c("1/5", "1/2", "1", "2", "5"))
ggheatmap_ten1000

pdf(file = "L1plot2.pdf", width = 9, height = 3)
ggarrange(ggheatmap_ten100, ggheatmap_ten300, ggheatmap_ten1000, ncol = 3, 
          align = "hv", widths = c(1, 1, 1), common.legend = TRUE, legend = "right")
dev.off()

#### densities

cnam <- c("5", "2", "1", "05", "02")
limits <- c(1, 2.4, 5, 10.5, 25, 1, 2.4, 5, 10.5, 25, 1, 2.4, 5, 10.5, 25)
resul <- matrix(0, ncol = 5, nrow = 4)
res_l <- list()
gre <- list()
df <- data.frame()

truegrid <- expand.grid(seq(-4, 4, length.out = 60), seq(-4, 4, length.out = 60))
trued <- cbind(truegrid, apply(truegrid, 1, dmvnorm))
names(trued) <- c("V1", "V2", "V3") 
truedl <- list()
mini <- 0
maxi <- 15
for(j in (mini + 1):maxi){
  truedl[[j - mini]] <- trued[abs(trued[,1]) < limits[j - mini] & abs(trued[,2]) < limits[j - mini],]
}

red <- '#FF0000FF'
grey <- '#303030'
blue <- '#3333B2'

for(i in (mini + 1):maxi){
  tdata2 <- tdata <- as.data.frame(using_list[[i]]$data)
  tcont <- as.data.frame(cbind(using_list[[i]]$grid, result[[i]][[1]]))
  
  gre[[i-mini]] <- ggscatter(data = tdata, x = "V1", y = "V2", color = grey, size = 3, alpha = 0.5) +
    stat_contour(geom="polygon", data = trued, mapping = aes(x = V1, y = V2, z = V3), alpha = 0.05, fill = "blue") +
    stat_contour(data = tcont, mapping = aes(x = V1, y = V2, z = V3), col = red, bins = 10) + theme_bw() +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) + border() + 
    coord_cartesian(xlim = c( -limits[i-mini],limits[i-mini]), ylim = c(-limits[i-mini],limits[i-mini])) 
}
pdf("DensShot_grid.pdf", width = 18, height = 10)
ggarrange(plotlist = gre, ncol = 5, nrow = 3)
dev.off()
