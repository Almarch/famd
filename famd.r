
library(FactoMineR)
library(agricolae)
library(ade4)
data(yacon)

yacon$height = cut(yacon$height,c(0,80,160,Inf))
n = nrow(yacon)

quanti = yacon[,c(7:10)]
quali = yacon[,c(1,3,6)]
p = ncol(quali)

# PCA ####
X_quanti = as.matrix(quanti - rep(colMeans(quanti), each = n))
s  <- sqrt(colSums(X_quanti^2) / n)
X_quanti <- X_quanti / rep(s, each = n) 

V = t(X_quanti) %*% X_quanti / n

eig <- eigen(V)
a <- eig$vectors
l <- eig$values

## normalization ####
for(j in 1:ncol(a)) a[,j] <- a[,j]  * sqrt(l[j])

## correlations ####
par(mfrow=c(1,2),mai=c(2,2,1,1)/2.54)
(lambda <- round(l*100/sum(l),1))
cercle=seq(0,2*pi, length=1000)

plot(c(-1,1),c(-1,1),type="n",
     xlab=paste0("Axe 1 : ",lambda[1],"% inertia"),
     ylab=paste0("Axe 2 : ",lambda[2],"% inertia"))
abline(v=0,col="black"); abline(h=0,col="black")
points(cos(cercle),sin(cercle),type="l",col="red")
arrows(0,0,a[,1], a[,2],angle=15,length=0.15)
text(a[,1]+0.05,a[,2]+.05,labels=colnames(quanti))

## individuals ####
coord = X_quanti %*% a
plot(coord[,1], coord[,2],
     xlab=paste0("Component 1 : ",lambda[1],"% inertia"),
     ylab=paste0("Component 2 : ",lambda[2],"% inertia"))
abline(v=0,col="black"); abline(h=0,col="black")

## FactoMineR ####
res = FactoMineR::PCA(yacon[,c(7:10)])

# MCA ####

X_quali = as.matrix(acm.disjonctif(quali))
D = diag(apply(X_quali,2,sum))

B = t(X_quali) %*% X_quali
eig = eigen(1/p * solve(D) %*% B) # Saporta 3rd edition p. 223

## remove the first eigen value (=1)
eig = list(vectors = eig$vectors[,-1],
           values = eig$values[-1])

## remove close-to_0 eigen values
is_not_0 = which(round(eig$values,10) != 0)
eig = list(vectors = eig$vectors[,is_not_0],
           values = eig$values[is_not_0])

a <- eig$vectors
l <- eig$values
(lambda <- round(l*100/sum(l),1))

## normalization ####
coord = c()
for(j in 1:ncol(a)){
  conv_norm = (1/(n*p) * t(a[,j]) %*% D %*% a[,j])[1,1]
  print(conv_norm)
  a[,j] = a[,j] * sqrt(l[j] / conv_norm)
  coord = cbind(
    coord,
    1/p * 1/sqrt(l[j]) * X_quali %*% a[,j]
  )
}

## individuals ####
par(mfrow=c(1,2),mai=c(2,2,1,1)/2.54)
plot(coord[,1], coord[,2],
     xlab=paste0("Component 1 : ",lambda[1],"% inertia"),
     ylab=paste0("Component 2 : ",lambda[2],"% inertia"))
abline(v=0,col="black"); abline(h=0,col="black")

## variables ####
plot(a[,1], a[,2],
     col = "red",pch = 4,cex=2,
     xlab=paste0("Component 1 : ",lambda[1],"% inertia"),
     ylab=paste0("Component 2 : ",lambda[2],"% inertia"))
abline(v=0,col="black"); abline(h=0,col="black")
text(a[,1]+0.05,a[,2]+.05,labels=colnames(X_quali), col = "red")

## FactoMineR ####
res = FactoMineR::MCA(yacon[,c(1,3,6)])

# FAMD ####

## Prepa data ####
X_quanti = as.matrix(quanti - rep(colMeans(quanti), each = n))
s  <- sqrt(colSums(X_quanti^2) / n)
X_quanti <- X_quanti / rep(s, each = n) 

X_quali = as.matrix(acm.disjonctif(quali))

D = diag(apply(X_quali,2,sum))

for(i in 1:ncol(X_quali)){
 proba = sqrt(sum(X_quali[,i]) / n)
 X_quali[,i] = X_quali[,i] / proba
}

#D = diag(apply(X_quali,2,sum))

X = cbind(X_quanti, X_quali)

## Eigen
V = t(X) %*% X / n

eig <- eigen(V)

## remove close-to_p eigen values
is_not_p = which(round(eig$values,10) != p)
eig = list(vectors = eig$vectors[,is_not_p],
           values = eig$values[is_not_p])

## remove close-to_0 eigen values
is_not_0 = which(round(eig$values,10) != 0)
eig = list(vectors = eig$vectors[,is_not_0],
           values = eig$values[is_not_0])

a <- eig$vectors
l <- eig$values
(lambda <- round(l*100/sum(l),1))

## normalization ####
prop = diag(D) / n
c = a[ncol(quanti) + 1:ncol(X_quali),]
coord = c()
for(j in 1:ncol(a)) {
  a[,j] <- a[,j] * sqrt(l[j])
  c[,j] <- a[ncol(quanti) + 1:ncol(X_quali), j] / sqrt(prop) * sqrt(l[j])
  coord = cbind(
    coord,
    1/sqrt(l[j]) * X %*% a[,j])
}

r = a[1:ncol(quanti),]

## correlations ####
par(mfrow=c(1,2),mai=c(2,2,1,1)/2.54)
(lambda <- round(l*100/sum(l),1))
cercle=seq(0,2*pi, length=1000)

plot(c(-1,1),c(-1,1),type="n",
     xlab=paste0("Axe 1 : ",lambda[1],"% inertia"),
     ylab=paste0("Axe 2 : ",lambda[2],"% inertia"))
abline(v=0,col="black"); abline(h=0,col="black")
points(cos(cercle),sin(cercle),type="l",col="red")
arrows(0,0,r[,1], r[,2],angle=15,length=0.15)
text(r[,1]+0.05,r[,2]+.05,labels=colnames(quanti))

## individuals + modalities
plot(coord[,1], coord[,2],
     xlab=paste0("Component 1 : ",lambda[1],"% inertia"),
     ylab=paste0("Component 2 : ",lambda[2],"% inertia"))
points(c[,1], c[,2],
       col = "red",pch = 4,cex=2,
       xlab=paste0("Axis 1 : ",lambda[1],"% inertia"),
       ylab=paste0("Axis 2 : ",lambda[2],"% inertia"))
abline(v=0,col="black"); abline(h=0,col="black")
text(c[,1]+0.05,c[,2]+.05,labels=colnames(X_quali), col = "red")

# FactoMineR
res = FactoMineR::FAMD(yacon[,c(1,3,6,7:10)])


