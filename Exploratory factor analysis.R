#Load Packages
library("corpcor"); 
library("GPArotation"); 
library("psych");
library(ellipse); 
library(nFactors); 
library(parallel); 
library(mvsf);

#Data load
mydata = read.csv("data.csv", header=FALSE,sep=";");
mydata.x <- data.frame(cbind(mydata$V5,mydata$V6,mydata$V7,mydata$V8,mydata$V9,mydata$V10,mydata$V11,mydata$V12,mydata$V13,mydata$V14,mydata$V15,mydata$V16,mydata$V17,
                             mydata$V18,mydata$V19,mydata$V20,mydata$V21,mydata$V22,mydata$V23,mydata$V24,mydata$V25,mydata$V26,mydata$V27,mydata$V28,mydata$V29,mydata$V30,mydata$V31,
                             mydata$V32,mydata$V33,mydata$V34,mydata$V35,mydata$V36,mydata$V37));
colnames(mydata.x) <- c("V5","V6","V7","V8","V9","V10","V11","V12","V13","V14","V15","V16","V17","V18","V19","V20","V21","V22","V23","V24","V25"
                        ,"V26","V27","V28","V29","V30","V31","V32","V33","V34","V35","V36","V37");

#Data prepare
mydata.x <- mydata.x[,c(1:4,7:22,29:30,33)]
colnames(mydata.x)

#Normality tests
mshapiro.test(t(mydata.x[1:179,]));
mvsf(t(mydata.x[1:179,]));

#KMO test
KMO(mydata.x[1:179,])

# Correlation matrix
matrix <- round(cor(mydata.x[1:179,]), digits=2); matrix;

#Numerical print
plotcorr(matrix, numbers = TRUE);

#Grephical represenation of the correlation matrix
colors <- c("#A50F15","#DE2D26","#FB6A4A","#FCAE91","#FEE5D9","white","#EFF3FF","#BDD7E7","#6BAED6","#3182BD","#08519C");
plotcorr(matrix, col=colors[5*matrix + 6]);

#Correlation matric only showing similarities above 0.35
colors <- c(
  "#A50F15","#A50F15","#A50F15","#A50F15","#DE2D26","#DE2D26","#DE2D26","#DE2D26","#FB6A4A","#FB6A4A",
  "#FB6A4A","#FB6A4A","#FB6A4A","white","white","white","white","white","white","white",
  "white",  
  "white","white","white","white","white","white","white","#6BAED6","#6BAED6","#6BAED6",
  "#6BAED6","#6BAED6","#3182BD","#3182BD","#3182BD","#3182BD","#08519C","#08519C","#08519C","#08519C");
plotcorr(matrix, col=colors[20*matrix + 21]);

#Remove variables with small correlation
mydata.y <- mydata.x[,c(8:13,15:22,23)]
colnames(mydata.y)
matrix <- round(cor(mydata.y[1:179,]), digits=2); matrix;

#Test for null values
cortest.bartlett(matrix,n=179);
#Test for multicolinearity must be larger than 0.00001 then there is no problem of multicolinearity
det(matrix, n = 179);

#Taking the final data and recheck KMO
mydata.final <- mydata.y[1:179,];
KMO(mydata.final);

#Eigenvalues Kaiser's criterion 
eigen(cor(mydata.final))

# Variance cummolative
fit <- princomp(mydata.final, cor=TRUE);
summary(fit) # print variance accounted for
# pc loadings
loadings(fit) 

#scree plot
plot(fit,type="lines") # scree plot
fit$scores # the principal components
biplot(fit); #graphical analysis

#Parallel analysis
ev <- eigen(cor(mydata.final)) # get eigenvalues
ap <- parallel(subject=nrow(mydata.final),var=ncol(mydata.final),rep=100,cent=.05);
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea);
plotnScree(nS);

#PCA
pca <- principal(mydata.final, nfactors=4, rotate="none");print.psych(pca, cut = 0.45, sort = TRUE);
pca <- principal(mydata.final, nfactors=4, rotate="varimax");print.psych(pca, cut = 0.45, sort = TRUE);
pca <- principal(mydata.final, nfactors=4, rotate="promax");print.psych(pca, cut = 0.45, sort = TRUE);

#ML - method 1
ml <- fa(mydata.final, nfactors=4, rotate="none", fm="ml");print.psych(ml, cut = 0.45, sort = TRUE);
ml <- fa(mydata.final, nfactors=4, rotate="varimax", fm="ml");print.psych(ml, cut = 0.45, sort = TRUE);
ml <- fa(mydata.final, nfactors=4, rotate="promax", fm="ml");print.psych(ml, cut = 0.45, sort = TRUE);
#ML - method 2
ml <- factanal(mydata.final, 4, rotation="none");print(ml, digits=2, cutoff=.45, sort=TRUE);
ml <- factanal(mydata.final, 4, rotation="varimax");print(ml, digits=2, cutoff=.45, sort=TRUE);
ml <- factanal(mydata.final, 4, rotation="promax");print(ml, digits=2, cutoff=.45, sort=TRUE);

#PAF
paf <- fa(mydata.final, nfactors=4, rotate="none", fm="pa");print.psych(paf, cut = 0.45, sort = TRUE);
paf <- fa(mydata.final, nfactors=4, rotate="varimax", fm="pa");print.psych(paf, cut = 0.45, sort = TRUE);
paf <- fa(mydata.final, nfactors=4, rotate="promax", fm="pa");print.psych(paf, cut = 0.45, sort = TRUE);

#Residuals check
residuals <- factor.residuals(matrix, paf$loadings)
residuals <- as.matrix(residuals[upper.tri(residuals)])
large.resid <- abs(residuals) > 0.05
sum(large.resid)
sum(large.resid)/nrow(residuals)
sqrt(mean(residuals^2))
hist(residuals)