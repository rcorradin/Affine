library(AFFINEpack)
library(doParallel)
library(reshape)
library(viridis)
library(ggplot2)
library(ggpubr)

# Get upper triangle of the matrix ----------------------------------------

get_upper_tri <- function(mat) {
  mat[lower.tri(mat)] <- NA
  return(mat)
}

# number of reps ----------------------------------------------------------

nrep <- 100

# t-STUD ------------------------------------------------------------------
# 100 ---------------------------------------------------------------------

data100 <- list()
for(i in 1:nrep){
  data100[[i]] <- rt(n = 100, df = 2)
}
grid <- seq(-100, 100, by = 0.1)

no_cores <- 6
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores)

a100 <- Sys.time()
result100 <- foreach(i = 1:nrep, .packages = c("AFFINEpack")) %dopar% {
  mod1 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
             data = 0.2 * data100[[i]], grid = 0.2 * grid, 
             m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
             m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)
  mod2 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 0.5 * data100[[i]], grid = 0.5 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)
  mod3 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 1 * data100[[i]], grid = 1 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)
  mod4 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 2 * data100[[i]], grid = 2 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)
  mod5 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 5 * data100[[i]], grid = 5 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)
  
  list(mod1, mod2, mod3, mod4, mod5)
}
b100 <- Sys.time()
b100 - a100

# t-STUD ------------------------------------------------------------------
# 300 ---------------------------------------------------------------------

data300 <- list()
for(i in 1:nrep){
  data300[[i]] <- rt(n = 300, df = 2)
}
grid <- seq(-100, 100, by = 0.1)

no_cores <- 6
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores)

a300 <- Sys.time()
result300 <- foreach(i = 1:nrep, .packages = c("AFFINEpack")) %dopar% {
  mod1 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 0.2 * data300[[i]], grid = 0.2 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)
  mod2 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 0.5 * data300[[i]], grid = 0.5 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)
  mod3 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 1 * data300[[i]], grid = 1 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)
  mod4 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 2 * data300[[i]], grid = 2 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)
  mod5 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 5 * data300[[i]], grid = 5 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)
  
  list(mod1, mod2, mod3, mod4, mod5)
}
b300 <- Sys.time()
b300 - a300

# t-STUD ------------------------------------------------------------------
# 1000 --------------------------------------------------------------------

data1000 <- list()
for(i in 1:nrep){
  data1000[[i]] <- rt(n = 1000, df = 2)
}
grid <- seq(-100, 100, by = 0.1)

no_cores <- 6
registerDoParallel(cores=no_cores)  
cl <- makeCluster(no_cores)

a1000 <- Sys.time()
result1000 <- foreach(i = 1:nrep, .packages = c("AFFINEpack")) %dopar% {
  mod1 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 0.2 * data1000[[i]], grid = 0.2 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)$dens
  mod2 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 0.5 * data1000[[i]], grid = 0.5 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)$dens
  mod3 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 1 * data1000[[i]], grid = 1 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)$dens
  mod4 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 2 * data1000[[i]], grid = 2 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)$dens
  mod5 <- MPU(niter = 7500, nburn = 2500, napprox = 100, 
              data = 5 * data1000[[i]], grid = 5 * grid, 
              m0 = 0, s0 = 1, a0 = 2, b0 = 1, 
              m1 = 0, k1 = 1, a1 = 2, b1 = 1, mass = 1, nupd = 250)$dens
  
  list(mod1, mod2, mod3, mod4, mod5)
}
b1000 <- Sys.time()
b1000 - a1000
save.image("IMAGE1000.RData")

# CLOSE -------------------------------------------------------------------

stopCluster(cl)

# PLOTS -------------------------------------------------------------------

ten100  <- array(0, dim = c(5, 5, 100))
ten300  <- array(0, dim = c(5, 5, 100))
ten1000 <- array(0, dim = c(5, 5, 100))
NGR   <- array(0, dim = c(5, 3, 100))
NGRVI <- array(0, dim = c(5, 3, 100))
cost <- c(0.2, 0.5, 1, 2, 5)

for(i in 1:100){
  for(j in 1:5){
    for(k in 1:5){
      ten100[j,k,i] <- sum(abs(result100[[i]][[j]] * cost[j] - 
                                 result100[[i]][[k]] * cost[k])) * (grid[3] - grid[2]) 
      ten300[j,k,i] <- sum(abs(result300[[i]][[j]] * cost[j] - 
                                 result300[[i]][[k]] * cost[k])) * (grid[3] - grid[2]) 
      ten1000[j,k,i] <- sum(abs(result1000[[i]][[j]] * cost[j] -
                                 result1000[[i]][[k]] * cost[k])) * (grid[3] - grid[2])
    }
    # NGR[j,1,i] <- mean(apply(result100[[i]][[j]]$clust, 1, function(x) length(unique(x))))
    # NGR[j,2,i] <- mean(apply(result300[[i]][[j]]$clust, 1, function(x) length(unique(x))))
    # NGR[j,3,i] <- mean(apply(result1000[[i]][[j]], 1, function(x) length(unique(x))))
  }
}
round(t(apply(NGR, c(1,2), mean)), digits = 2)

# Plots

upper_tri_ten100 <- get_upper_tri(apply(ten100, c(1,2), mean))/max(apply(ten100, c(1,2), mean))
melted_ten100 <- melt(upper_tri_ten100, na.rm = TRUE)
ggheatmap_ten100 <- ggplot(melted_ten100, aes(X2, X1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "D", name = "", direction=-1, na.value = "white", limits = c(0,1), guide = F) +
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

upper_tri_ten300 <- get_upper_tri(apply(ten300, c(1,2), mean))/max(apply(ten100, c(1,2), mean))
melted_ten300 <- melt(upper_tri_ten300, na.rm = TRUE)
ggheatmap_ten300 <- ggplot(melted_ten300, aes(X2, X1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "D", name = "", direction=-1, na.value = "white", limits = c(0,1), guide = FALSE) +
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

upper_tri_ten1000 <- get_upper_tri(apply(ten1000, c(1,2), mean))/max(apply(ten100, c(1,2), mean))
melted_ten1000 <- melt(upper_tri_ten1000, na.rm = TRUE)
ggheatmap_ten1000 <- ggplot(melted_ten1000, aes(X2, X1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "D", name = "", direction=-1, na.value = "white", limits = c(0,1)) +
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

pdf(file = "L1plot2_tsutdent.pdf", width = 9, height = 3)
ggarrange(ggheatmap_ten100, ggheatmap_ten300, ggheatmap_ten1000, ncol = 3, 
          align = "hv", widths = c(1, 1, 1), common.legend = TRUE, legend = "right")
dev.off()

