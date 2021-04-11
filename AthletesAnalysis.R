# Athletes Analysis

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(GGally)
library(dplyr) 
library(tidyverse)
library(dendextend)
library(MASS)
library(lubridate)
library(Matrix)
library(ggpubr)
library(factoextra)
library(ggfortify)
library(ISLR)
library(pls)
library(corrplot)
library(CCA)
library(CCP)

# reading in data
ais0 = read.csv("ais.csv")
      # make Sex a factor
ais0 = ais0 %>%
   mutate( Sex = factor( Sex ))

# subset without Sex and LBM
ais.q1 = ais0[c(-9,-7)]

# parallel coordinate plot (all columns except group by column)
ggparcoord(ais0, columns=c(1,2,3,4,5,6,7,8,10,11,12),
        groupColumn = 9, scale="globalminmax", title = "Parallel Coordinate Plot")


### PCA ###
pca <- prcomp(ais.q1)
      # shows PC proportions/cumulative proportion
summary(pca)
      # shows variable contribution per PC
pca

# eigenvalue/scree plot
plot(pca$sdev^2, type="b", main="Eigenvalues")

 # cumulative variance plot
plot(cumsum(pca$sdev^2 / sum(pca$sdev^2)), type="b", main="Cumulative Variance")
      # alternative way - non-ratio
fviz_screeplot(pca, addlabels = TRUE)

# covariance matrix of data
cov_mat = cov(ais.q1)
      # eigenvalues
eigen(cov_mat)$value
      # largest eigenvalue
max(eigen(cov_mat)$value)

# PC score plots of raw data
autoplot(pca, x=1, y=2, data = ais0, colour = "Sex", main = "PC1 against PC2 - raw data")
autoplot(pca, x=2, y=3, data = ais0, colour = "Sex", main = "PC2 against PC3 - raw data")

# scale data
ais.q1.sc = scale(ais.q1, center = TRUE, scale = TRUE)

# PCA of scaled data
pca2 <- prcomp(ais.q1.sc)
summary(pca2)
pca2

# variable contributions in decreasing order
      # raw data - PC1 and PC2
sort(abs(pca$rotation[,1]), decreasing = TRUE)
sort(abs(pca$rotation[,2]), decreasing = TRUE)
      # scaled data - PC1 and PC2
sort(abs(pca2$rotation[,1]), decreasing = TRUE)
sort(abs(pca2$rotation[,2]), decreasing = TRUE)

# PC score plots of scaled data
autoplot(pca2, x=1, y=2, data = ais0, colour = "Sex", main = "PC1 against PC2 - scaled data")
autoplot(pca2, x=2, y=3, data = ais0, colour = "Sex", main = "PC2 against PC3 - scaled data")

# extra - density plots
      # raw data - PC1 and PC2
plot(density(pca$x[,1]), xlab = "PC1")
plot(density(pca$x[,2]), xlab = "PC2")
      # scaled data - PC1 and PC2
plot(density(pca2$x[,1]), xlab = "PC1")
plot(density(pca2$x[,2]), xlab = "PC2")


### PCR ###
# pca = prcomp(data=) #remove the response variable from data
pcaScores = as.data.frame(pca$x)
pcaScores$LBM = ais0$LBM #add response variable to the pcaScores dataframe
pcr.fwd.r = pcr.null.r = lm(LBM ~., data=pcaScores) #null model if needed
PCAnames = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10") # List the amount you need
modelsList = list(numeric(10))
for (i in 1:10){
   modelsList[[i]] = paste("LBM ~", paste(PCAnames[1:i],collapse='+'))
}
modelsList = lapply(modelsList, as.formula)
modelsResults = lapply(modelsList, function(x) lm(x, data = pcaScores))
modelsSummary = lapply(modelsResults, summary)
modelsSummary
residualSD = data.frame()
for (i in 1:10){
   temp.r = data.frame(Cluster_num=i, sigma=modelsSummary[[i]]$sigma)
   residualSD = rbind.data.frame(residualSD, temp.r)
}
residualSD
plot(residualSD, type = "b", main = "Residual Standard Deviations")


### CORRELATION / CCA ###
# subset without Sex
ais.q2 = ais0[c(-9)]

# two subsets
X = ais.phys = ais0[c(1,2,6,7,10,12)]
Y = ais.sero = ais0[c(3,4,5,8,11)]

# between variate correlations (between correlation matrix)
betcor = cor(X, Y)
betcor
      # strongest correlation
max(abs(betcor))
# plot of most strongly correlated variables
plot(X$PctBfat, Y$Hg, main = "Most Correlated Variables")

# also: non CCA between covariance matrix
betcov = cov(X, Y)
betcov

# CCA model
cca <- cancor(X, Y)
      # canonical correlations in order of strength
cca$cor

# extra: compute canonical loadings
cca.loadings <- comput(X, Y, cca)
      # plot U and V scores
plot(cca.loadings$corr.X.xscores[,1], type = "b", main = "U Scores")
plot(cca.loadings$corr.Y.yscores[,1], type = "b", main = "V Scores")

# canonical covariate pairs
CC1_X <- as.matrix(X) %*% cca$xcoef[, 1]
CC1_Y <- as.matrix(Y) %*% cca$ycoef[, 1]

CC2_X <- as.matrix(X) %*% cca$xcoef[, 2]
CC2_Y <- as.matrix(Y) %*% cca$ycoef[, 2]

CC3_X <- as.matrix(X) %*% cca$xcoef[, 3]
CC3_Y <- as.matrix(Y) %*% cca$ycoef[, 3]

CC4_X <- as.matrix(X) %*% cca$xcoef[, 4]
CC4_Y <- as.matrix(Y) %*% cca$ycoef[, 4]

CC5_X <- as.matrix(X) %*% cca$xcoef[, 5]
CC5_Y <- as.matrix(Y) %*% cca$ycoef[, 5]
      # plots of each pair of canonical correlation scores
plot(CC1_X, CC1_Y, colour = ais0$Sex)
plot(CC2_X, CC2_Y)
plot(CC3_X, CC3_Y)
plot(CC4_X, CC4_Y)
plot(CC5_X, CC5_Y)

# hypothesis testing
Tk = function( k, n, d1, d2, vv, sig ){
   Tkout = - ( n - ( d1 + d2 + 3 )/2 ) * log( prod ( 1 - vv[1 ]^2 ) )
   dof = ( d1 - k ) * ( d2 - k )
   pval = pchisq( Tkout, df = dof, lower.tail = FALSE )
   crit = qchisq(sig, df = dof, lower.tail=TRUE)
   list( Tkout = Tkout, pval = pval, dof = dof, crit = crit ) }
Tk( k = 1, n = 202, d1 = 6, d2 = 5, vv = cca$cor[2], sig = 0.02 )
Tk( k = 2, n = 202, d1 = 6, d2 = 5, vv = cca$cor[3], sig = 0.02 )
Tk( k = 3, n = 202, d1 = 6, d2 = 5, vv = cca$cor[4], sig = 0.02 )
Tk( k = 4, n = 202, d1 = 6, d2 = 5, vv = cca$cor[5], sig = 0.02 )
# above for k from 1-4, change k and vv number to k+1 each time


### K-MEANS CLUSTERING ###
# subset with variables in required order
ais.q3 = cbind(ais.phys, ais.sero)
      # scale data
ais.q3.sc = scale(ais.q3)
      # split data
ais.phys.sc = ais.q3.sc[,1:6]
ais.sero.sc = ais.q3.sc[,7:11]

## for ais.phys.sc
# all kmeans in a list
kmeans.list = list()
for(i in 1:4){ # <- choose number of k here
   kmeans.list[[i]] = kmeans(ais.phys.sc, centers=i, nstart= 25)
}
# kmeans.list

# membership table - all k
clust.mem = data.frame()
for(i in 1:4){
   temp.info = data.frame(Total_Clusters=i, Cluster_Label=kmeans.list[[i]]$cluster)
   clust.mem = rbind.data.frame(clust.mem, temp.info)
}
xtabs(~Cluster_Label+Total_Clusters, data=clust.mem)

# cluster alignment with factor per k
table(ais0$Sex, kmeans.list[[2]]$cluster) 
table(ais0$Sex, kmeans.list[[3]]$cluster) 
table(ais0$Sex, kmeans.list[[4]]$cluster) 

## for ais.sero.sc
# all kmeans in a list
kmeans.list = list()
for(i in 1:4){ # <- choose number of k here
   kmeans.list[[i]] = kmeans(ais.sero.sc, centers=i, nstart= 25)
}
# kmeans.list

# membership table - all k
clust.mem = data.frame()
for(i in 1:4){
   temp.info = data.frame(Total_Clusters=i, Cluster_Label=kmeans.list[[i]]$cluster)
   clust.mem = rbind.data.frame(clust.mem, temp.info)
}
xtabs(~Cluster_Label+Total_Clusters, data=clust.mem)

# cluster alignment with factor per k
table(ais0$Sex, kmeans.list[[2]]$cluster) 
table(ais0$Sex, kmeans.list[[3]]$cluster) 
table(ais0$Sex, kmeans.list[[4]]$cluster) 

# extra - within and between cluster variabilities
within.df = data.frame()
between.df = data.frame()
for(i in 1:4){
   temp.w = data.frame(Total_Clusters=i, Within_Variance=kmeans.list[[i]]$tot.withinss)
   temp.b = data.frame(Total_Clusters=i, Between_Variance=kmeans.list[[i]]$betweenss)
   within.df = rbind.data.frame(within.df, temp.w)
   between.df = rbind.data.frame(between.df, temp.b)
}
within.df
plot(within.df, type = "b", main = "Within Cluster Variabiliy")
between.df
plot(between.df, type = "b", main = "Between Cluster Variabiliy")

# error rate:
# manual way for k means 2 cluster case -> (errors/total)*100
err = round(((9+15)/(9+15+85+93))*100,2)
err

## LDA ###
# with uniform priors - depends on number of factors (2 factors here)
# (scaling doesn't matter)
# for ais.phys 
lda_phys =  lda(Sex~. , data = ais0[c(1,6,7,9,10,12)], prior=c(1,1)/2)
      # confusion matrix
      # non-diagonals are misclassified:
      # column value is misclassified as row value
lda_phys_pred <- predict(lda_phys)$class
table(lda_phys_pred, ais0$Sex)
      # error rate
round(mean(lda_phys_pred != ais0$Sex)*100, 2)

# for ais.sero
lda_sero =  lda(Sex~. , data = ais0[c(1,4,5,8,9,11)], prior=c(1,1)/2)
      # confusion matrix
lda_sero_pred <- predict(lda_sero)$class
table(lda_sero_pred, ais0$Sex)
      # error rate
round(mean(lda_sero_pred != ais0$Sex)*100, 2)


### LOGISTIC REGRESSION ###
# based on glm
# (scaling doesn't matter)
## for ais.phys
log_phys = glm(Sex~. , data = ais0[c(1,6,7,9,10,12)], family = binomial)
log_phys_pred <- predict(log_phys, type = "link")
      # change m and f to class labels, default threshold of 0.5
phys_pred_class <- ifelse(log_phys_pred > 0.5, "m", "f")
      # confusion matrix - may have to swap "m" and "f" values above?
table(phys_pred_class, ais0$Sex)
      # error rate
round(mean(phys_pred_class != ais0$Sex)*100, 2)

## for ais.sero
log_sero = glm(Sex~. , data = ais0[c(1,4,5,8,9,11)], family = binomial)
log_sero_pred <- predict(log_sero, type = "link")
      # change m and f to class labels, default threshold of 0.5
sero_pred_class <- ifelse(log_sero_pred > 0.5, "m", "f")
      # confusion matrix - may have to swap "m" and "f" values above?
table(sero_pred_class, ais0$Sex)
      # error rate
round(mean(sero_pred_class != ais0$Sex)*100, 2)

# extra - probability histogram: do for each pred model
      # convert from log to standard
odds = exp(log_sero_pred)
      # convert from odds to probabilities
probabilities = odds/(odds+1)
      # histogram
hist(probabilities)
