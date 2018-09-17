library(AFFINEpack)
library(reshape)
library(viridis)
library(ggplot2)
library(ggpubr)
library(mcclust.ext)
library(clv)

#-----------------------#
#  Import data + clean  #
#-----------------------#

data0 <- read.table('dataXY2.dat', skip = 2)
head(data0)
plot(data0[,c(7,8,10,11)])

#----------
# coherently with the work done
# by Ibata et al. (2011)
# we choose for the stars with multiple 
# obs the one with the lowest 
# estimated velocity uncertainty

max_ind <- max(data0[,1])
dataclean <- data0[data0[,1] == 1,]
for(j in 2:max_ind){
  temp <- data0[data0[,1]==j,]
  if(length(dim(temp)) == 2){temp = temp[which.min(temp[,5]),]}
  dataclean <- rbind(dataclean, temp)
}
plot(dataclean[,c(4,8,10,11)])

# take only variables V, FeH and R
# discard other variables

dataNS <- dataclean[dataclean[,8] > -2.5, c(4,8,10,11)]
plot(dataNS)
data <- scale(dataNS)
range(data[,3])
plot(as.data.frame(data))

grid <- expand.grid(seq(min(data[,1]) - 0.2, max(data[,1]) + 0.2, length.out = 40),
                    seq(min(data[,2]) - 0.2, max(data[,2]) + 0.2, length.out = 40),
                    seq(min(data[,3]) - 0.2, max(data[,3]) + 0.2, length.out = 40),
                    seq(min(data[,4]) - 0.2, max(data[,4]) + 0.2, length.out = 40))

#--------------------#
#  Model Estimation  #
#--------------------#

set.seed(42)
m_0_19 <- DPmixMulti(nsim = 25000,
                     nburn = 5000, 
                     napprox = 100,
                     data = as.matrix(data), 
                     grid = as.matrix(grid), 
                     mu_start = c(0,0,0,0), 
                     Lambda_start = diag(1,4), 
                     theta = 1, 
                     m0 = c(0,0,0,0), 
                     B0 = diag(1,4), 
                     nu0 = 26, 
                     sigma = diag(21,4), 
                     b1 = 6, 
                     B1 = diag(15,4), 
                     m1 = c(0,0,0,0), 
                     k1 = 1,
                     t1 = 1, 
                     t2 = 5.26, 
                     nupd = 100, 
                     fix = FALSE,
                     dep = TRUE)

#-----------#
#   graph   #
#-----------#

grey <- '#868686FF'


data <- as.data.frame(dataNS)
names(data) <- c("V", "FeH", "X", "Y")
gridS <- cbind(grid[,1]*sqrt(var(data[,1])) + mean(data[,1]), 
               grid[,2]*sqrt(var(data[,2])) + mean(data[,2]), 
               grid[,3]*sqrt(var(data[,3])) + mean(data[,3]), 
               grid[,4]*sqrt(var(data[,4])) + mean(data[,4]))
temp <- as.data.frame(cbind(gridS, m_0_19[[1]]))
names(temp) <- c("V1", "V2", "V3", "V4", "V5")

##### NO DENSITY

spVFeHn <- ggscatter(data, x = "V", y = "FeH", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression(V)) + ylab(expression("[Fe/H]"))
spVXn <- ggscatter(data, x = "V", y = "X", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression(V)) + ylab(expression(Y[1]))
spVYn <- ggscatter(data, x = "V", y = "Y", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression(V)) + ylab(expression(Y[2]))
spFeHXn <- ggscatter(data, x = "FeH", y = "X", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression("[Fe/H]")) + ylab(expression(Y[1]))
spFeHYn <- ggscatter(data, x = "FeH", y = "Y", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression("[Fe/H]")) + ylab(expression(Y[2]))
spXYn <- ggscatter(data, x = "X", y = "Y", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression(Y[1])) + ylab(expression(Y[2]))

pdf("pic_joint.pdf", width = 10, height = 8)
ggarrange(spVFeHn,  NULL, NULL, 
          spVXn, spFeHXn, NULL, 
          spVYn, spFeHYn, spXYn,
          ncol = 3, nrow = 3,  align = "hv", 
          widths = c(1, 1, 1), heights = c(1, 1, 1))
dev.off()

##### DENSITY

tdata1 <- aggregate(temp$V5, by = list(temp$V1, temp$V2), FUN = sum)
tdata2 <- aggregate(temp$V5, by = list(temp$V1, temp$V3), FUN = sum)
tdata3 <- aggregate(temp$V5, by = list(temp$V1, temp$V4), FUN = sum)
tdata4 <- aggregate(temp$V5, by = list(temp$V2, temp$V3), FUN = sum)
tdata5 <- aggregate(temp$V5, by = list(temp$V2, temp$V4), FUN = sum)
tdata6 <- aggregate(temp$V5, by = list(temp$V3, temp$V4), FUN = sum)

clB <- "#FF0000FF"
black <- '#303030'
spVFeH <- ggscatter(data, x = "V", y = "FeH", color = black, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata1, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7, bins = 10)+ 
  xlab(expression(V)) + ylab(expression("[Fe/H]")) + theme_bw()
spVX <- ggscatter(data, x = "V", y = "X", color = black, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata2, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7, bins = 10) + 
  xlab(expression(V)) + ylab(expression(Y[1])) + theme_bw()
spVY <- ggscatter(data, x = "V", y = "Y", color = black, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata3, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7, bins = 10) + 
  xlab(expression(V)) + ylab(expression(Y[2])) + theme_bw()
spFeHX <- ggscatter(data, x = "FeH", y = "X", color = black, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata4, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7, bins = 10) + 
  xlab(expression("[Fe/H]")) + ylab(expression(Y[1])) + theme_bw()
spFeHY <- ggscatter(data, x = "FeH", y = "Y", color = black, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata5, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7, bins = 10) + 
  xlab(expression("[Fe/H]")) + ylab(expression(Y[2])) + theme_bw()
spXY <- ggscatter(data, x = "X", y = "Y", color = black, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata6, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7, bins = 10) + 
  xlab(expression(Y[1])) + ylab(expression(Y[2])) + theme_bw()

pdf("pic_cont_m_grid.pdf", width = 10, height = 8)
ggarrange(spVFeH, NULL, NULL, 
          spVX, spFeHX, NULL, 
          spVY, spFeHY, spXY,
          ncol = 3, nrow = 3,  align = "hv", 
          widths = c(1, 1, 1), heights = c(1, 1, 1))
dev.off()

#----------------#
# Best Partition #
#----------------#

psm=comp.psm(m_0_19[[2]] + 1)
data.VI=minVI(psm,m_0_19[[2]] + 1,method=("all"),include.greedy=TRUE)
summary(data.VI)

#### Best VI and color

datac <- cbind(data, data.VI$cl[1,], as.factor(ifelse(data.VI$cl[1,] == 1, 1, 2)))
datac[,5] <- as.factor(datac[,5])
names(datac) <- c("V", "FeH" ,"X", "Y", "group", "shap")
cbPalette <- c("#303030", "#ff9d00", "#FF0000FF", "#006813", "#3333B2")

spVFeH <- ggscatter(datac, x = "V", y = "FeH", shape = "shap", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  border() + xlab(expression(V)) + ylab(expression("[Fe/H]")) + theme_bw() + theme(legend.position="none") 
spFeHX <- ggscatter(datac, x = "FeH", y = "X", shape = "shap", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  border() + xlab(expression("[Fe/H]")) + ylab(expression(Y[1])) + theme_bw() + theme(legend.position="none") 
spFeHY <- ggscatter(datac, x = "FeH", y = "Y", shape = "shap", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  border() + xlab(expression("[Fe/H]")) + ylab(expression(Y[2])) + theme_bw() + theme(legend.position="none") 
spVX <- ggscatter(datac, x = "V", y = "X", shape = "shap", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  border() + xlab(expression(V)) + ylab(expression(Y[1])) + theme_bw() + theme(legend.position="none") 
spVY <- ggscatter(datac, x = "V", y = "Y", shape = "shap", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  border() + xlab(expression(V)) + ylab(expression(Y[2])) + theme_bw() + theme(legend.position="none") 
spXY <- ggscatter(datac, x = "X", y = "Y", shape = "shap", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  border() +xlab(expression(Y[1])) + ylab(expression(Y[2])) + theme_bw() + theme(legend.position="none") 

pdf("pic_W_best_m_grid.pdf", width = 10, height = 8)
ggarrange(spVFeH, NULL, NULL, 
          spVX, spFeHX, NULL, 
          spVY, spFeHY, spXY,
          ncol = 3, nrow = 3,  align = "hv", 
          widths = c(1, 1, 1), heights = c(1, 1, 1))
dev.off()

#######################################################################

mat <- comp.psm(m_0_19[[2]] + 1)
reorder_mat <- function(mat) {
  dd <- as.dist((1 - mat) / 2)
  hc <- hclust(dd)
  mat <- mat[hc$order, hc$order]
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(mat) {
  mat[lower.tri(mat)] <- NA
  return(mat)
}
mat <- reorder_mat(mat)
upper_tri <- get_upper_tri(mat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggheatmap <- ggplot(melted_cormat, aes(X2, X1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "viridis", name = "Posterior\nsimilarity", direction=-1, na.value = "white") +
  theme_minimal() + 
  coord_fixed() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_colorbar(
    barwidth = 7,
    barheight = 1,
    title.position = "top",
    title.hjust = 0.5
  ))

pdf("heat-m.pdf", width = 8, height = 8)
ggheatmap
dev.off()
