library(BNPmix)
library(reshape)
library(viridis)
library(ggplot2)
library(ggpubr)
library(mcclust.ext)

#-----------------------#
#  Import data + clean  #
#-----------------------#

data0 <- read.table('dataXY.dat', skip = 2)
head(data0)

# take only variables V, FeH and R
# discard other variables

dataNS <- data0[data0[,8] > -2.5, c(7,8,10,11)]
nrow(data)
plot(data)
data <- scale(dataNS)
range(data[,3])

grid <- expand.grid(seq(min(data[,1]) - 0.2, max(data[,1]) + 0.2, length.out = 40),
                    seq(min(data[,2]) - 0.2, max(data[,2]) + 0.2, length.out = 40),
                    seq(min(data[,3]) - 0.2, max(data[,3]) + 0.2, length.out = 40),
                    seq(min(data[,4]) - 0.2, max(data[,4]) + 0.2, length.out = 40))

#--------------------#
#  Model Estimation  #
#--------------------#

set.seed(42)
m_2_42 <- DPmixMulti(nsim = 25000,
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
                    B1 = diag(10,4), 
                    m1 = c(0,0,0,0), 
                    k1 = 1,
                    t1 = 1, 
                    t2 = 2.41745, 
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
temp <- as.data.frame(cbind(gridS, m_2_42[[1]]))
names(temp) <- c("V1", "V2", "V3", "V4", "V5")

##### NO DENSITY

spVFeHn <- ggscatter(data, x = "V", y = "FeH", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression(V)) + ylab(expression(Fe/H))
spVXn <- ggscatter(data, x = "V", y = "X", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression(V)) + ylab(expression(Y[1]))
spVYn <- ggscatter(data, x = "V", y = "Y", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression(V)) + ylab(expression(Y[2]))
spFeHXn <- ggscatter(data, x = "FeH", y = "X", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression(Fe/H)) + ylab(expression(Y[1]))
spFeHYn <- ggscatter(data, x = "FeH", y = "Y", color = grey, size = 3, alpha = 0.5)+ border() + 
  xlab(expression(Fe/H)) + ylab(expression(Y[2]))
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

clB <- "#8F0000" 
spVFeH <- ggscatter(data, x = "V", y = "FeH", color = grey, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata1, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7)+ 
  xlab(expression(V)) + ylab(expression(Fe/H))
spVX <- ggscatter(data, x = "V", y = "X", color = grey, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata2, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7) + 
  xlab(expression(V)) + ylab(expression(Y[1]))
spVY <- ggscatter(data, x = "V", y = "Y", color = grey, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata3, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7) + 
  xlab(expression(V)) + ylab(expression(Y[2]))
spFeHX <- ggscatter(data, x = "FeH", y = "X", color = grey, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata4, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7) + 
  xlab(expression(Fe/H)) + ylab(expression(Y[1]))
spFeHY <- ggscatter(data, x = "FeH", y = "Y", color = grey, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata5, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7) + 
  xlab(expression(Fe/H)) + ylab(expression(Y[2]))
spXY <- ggscatter(data, x = "X", y = "Y", color = grey, size = 3, alpha = 0.5)+ border() +
  stat_contour(data = tdata6, mapping = aes(x = Group.1, y = Group.2, z = x), col = clB, alpha = 0.7) + 
  xlab(expression(Y[1])) + ylab(expression(Y[2]))

pdf("pic_cont_m.pdf", width = 10, height = 8)
ggarrange(spVFeH, NULL, NULL, 
          spVX, spFeHX, NULL, 
          spVY, spFeHY, spXY,
          ncol = 3, nrow = 3,  align = "hv", 
          widths = c(1, 1, 1), heights = c(1, 1, 1))
dev.off()

#----------------#
# Best Partition #
#----------------#

psm=comp.psm(m_2_42[[2]] + 1)
data.VI=minVI(psm,m_2_42[[2]] + 1,method=("all"),include.greedy=TRUE)
summary(data.VI)

#### Best VI and color

datac <- cbind(data, data.VI$cl[1,])
datac[,5] <- as.factor(datac[,5])
names(datac) <- c("V", "FeH" ,"X", "Y", "group")
cbPalette <- c("#999999", "#3333B2", "#8F0000")

spVFeH <- ggscatter(datac, x = "V", y = "FeH",  shape = "group", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  theme(legend.position="none") + border() + 
  xlab(expression(V)) + ylab(expression(Fe/H))
spFeHX <- ggscatter(datac, x = "FeH", y = "X", shape = "group", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  theme(legend.position="none") + border() + 
  xlab(expression(Fe/H)) + ylab(expression(Y[1]))
spFeHY <- ggscatter(datac, x = "FeH", y = "Y", shape = "group", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  theme(legend.position="none") + border() + 
  xlab(expression(Fe/H)) + ylab(expression(Y[2]))
spVX <- ggscatter(datac, x = "V", y = "X", shape = "group", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  theme(legend.position="none") + border() + 
  xlab(expression(V)) + ylab(expression(Y[1]))
spVY <- ggscatter(datac, x = "V", y = "Y", shape = "group", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  theme(legend.position="none") + border() + 
  xlab(expression(V)) + ylab(expression(Y[2]))
spXY <- ggscatter(datac, x = "X", y = "Y", shape = "group", color = "group", palette = cbPalette, size = 3, alpha = 0.5) + 
  theme(legend.position="none") + border() +
  xlab(expression(Y[1])) + ylab(expression(Y[2]))
pdf("pic_W_best_m.pdf", width = 10, height = 8)
ggarrange(spVFeH, NULL, NULL, 
          spVX, spFeHX, NULL, 
          spVY, spFeHY, spXY,
          ncol = 3, nrow = 3,  align = "hv", 
          widths = c(1, 1, 1), heights = c(1, 1, 1))
dev.off()

#-----------------#
#     Heatmap     #
#-----------------#

reorder_mat <- function(mat) {
  # Use correlation between variables as distance
  dd <- as.dist((1 - mat) / 2)
  hc <- hclust(dd)
  mat <- mat[hc$order, hc$order]
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(mat) {
  mat[lower.tri(mat)] <- NA
  return(mat)
}

psm_ord <- reorder_mat(comp.psm(m_2_42[[2]] + 1))
upper_tri <- get_upper_tri(psm_ord)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggheatmap <- ggplot(melted_cormat, aes(X2, X1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "viridis", name = "Posterior\nsimilarity", direction=-1, na.value = "white") +
  theme_minimal() + # minimal theme
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
