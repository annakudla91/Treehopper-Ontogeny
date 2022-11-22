### Repeatability Analysis ###
#Part I - 10% (36 specimens) of the dataset was landmarked a second time
#Part II - 5% (18 of the 36 specimens selected for Part I) of the dataset was 
#photographed a second time and landmarked

## Part I
setwd("~/Dropbox (Duke Bio_Ea)/Projects/Chapter 2 - Treehopper Ontogeny/R Analyses/All Instar Stages/Repeatability Analysis")
library(geomorph)
library(Morpho)
library(shapes)
library(ggplot2)
library(GGally)
library(lmerTest)
library(lme4)
library(emmeans)
library(multcomp)
library(lattice)
library(Momocs)

## Load in the data
# Data For Part I #
Original.36.land<-scan("OriginalLandmarksRepeatability-36-final.reorder.TPS", what = "character")
Repeat.36.land<-scan("Repeatability-S3.tps", what = "character")

# Data For Part II #
Photo1.18<-scan("Repeatability-Rep2_S3.TPS", what = "character")
Photo2.18<-scan("Repeatability-Rep3_S3.TPS", what = "character")

# Preparing Landmarks For Part I #
# 52 elements for each specimen in the dataset and 36 specimens
1872/52
# identify which lines of the vector belong to scale
#Using regular sequence, we can get image names and coordinates.
Scale.Set1<-Original.36.land[c(1:36)*52]
Scale.Set1
scale.Set1<-sub("SCALE=", "", Scale.Set1)
scale.Set1
mode(scale.Set1)<-"numeric"
is.numeric(scale.Set1)
Original.36.land[1:150]

# generate an empty array to store the 36 configurations (specimens) for Rep 1
n<-36
Set1.land<-array(NA, dim=c(24,2,36)) 
for (i in 1:n){Set1.land[,,i]<-matrix(Original.36.land[(i*52-50):(i*52-3)],24,2,byrow=T)}

Set1.land[,,12]
mode(Set1.land)<-"numeric"

# set the scale
Set1.landmarks<-Set1.land
for (i in 1:36)
{Set1.landmarks[,,i]<-(Set1.land[,,i]*(scale.Set1[i]))}
Set1.landmarks<-Set1.landmarks[1:24,,]
Set1.landmarks

dim(Set1.landmarks)
summary(Set1.landmarks)

Scale.Set2<-Repeat.36.land[c(1:36)*52]
Scale.Set2
scale.Set2<-sub("SCALE=", "", Scale.Set2)
scale.Set2
mode(scale.Set2)<-"numeric"
is.numeric(scale.Set2)

# generate an empty array to store the 36 configurations (specimens) for Rep 2
n<-36
Set2.land<-array(NA, dim=c(24,2,36)) 
for (i in 1:n){Set2.land[,,i]<-matrix(Repeat.36.land[(i*52-50):(i*52-3)],24,2,byrow=T)}

Set2.land[,,12]
mode(Set2.land)<-"numeric"

# want to use this loop to set the scale, scale is 1mm=0.001 pixels or so
# set the scale
Set2.landmarks<-Set2.land
for (i in 1:36)
{Set2.landmarks[,,i]<-(Set2.land[,,i]*(scale.Set2[i]))}
Set2.landmarks<-Set2.landmarks[1:24,,]
summary(Set2.landmarks)

dim(Set2.landmarks)
plot(Set1.landmarks)
plot(Set2.landmarks)
Set2.landmarks[,,14]
Set1.landmarks[,,14]
plot(Set1.landmarks[,,36])

p<-72
RepsAll<-array(NA, dim=c(24,2,p))
RepsAll[,,1:36]<-Set1.landmarks
RepsAll[,,37:72]<-Set2.landmarks
dim(RepsAll)
plot(RepsAll)

# Procrustes Fit 
Rep.sup<-gpagen(RepsAll)
plot(Rep.sup)

# Reps of individual from the same photo, landmarked 2x
ind1<-1:36
ind2<-1:36
inds<-c(ind1, ind2)
inds<-as.factor(inds)
size<-Rep.sup$Csize

# model 1, size by individuals
mod1<-summary(aov(size~inds))
s2within<- mod1[[1]][2,3]
s2among<- (mod1[[1]][1,3] - s2within)/2
#We compute the % of measurement error for size.
s2within/(s2within+s2among) # 0.46% error in size by landmarking

shape<-Rep.sup$coords
S<-two.d.array(shape)
#For shape, things are similar, we must estimate SS and crossproduct matrices first.
# model 2, shape by individuals 
mod2<-lm(S~inds)
SSres<-crossprod(residuals(mod2)) 
SSeffect<-crossprod(scale(mod2$fitted, scale=F))
tracewithin<-sum(diag(SSres))/ mod2$df.residual
traceamong<-(sum(diag(SSeffect))/(length(unique(inds))-1) - tracewithin)/2
#We compute the % of shape measurement error.
tracewithin/(tracewithin+traceamong) # 5% shape measurement error


## Part II - Different Images
subset<-c(1:36)
sample(subset, size = 18)

# Set 1 Rep 2
Scale.1.rep2<-Photo1.18[c(1:18)*52]
Scale.1.rep2
scale.1.rep2<-sub("SCALE=", "", Scale.1.rep2)
scale.1.rep2
mode(scale.1.rep2)<-"numeric"
is.numeric(scale.1.rep2)

# generate an empty array to store the 18 configurations (specimens) for Rep 2 of Test II
n<-18
Set1.2.land<-array(NA, dim=c(24,2,18)) 
for (i in 1:n){Set1.2.land[,,i]<-matrix(Photo1.18[(i*52-50):(i*52-3)],24,2,byrow=T)}

Set1.2.land[,,12]
mode(Set1.2.land)<-"numeric"

# set the scale
Set1.2.landmarks<-Set1.2.land
for (i in 1:18)
{Set1.2.landmarks[,,i]<-(Set1.2.land[,,i]*(scale.1.rep2[i]))}
Set1.2.landmarks<-Set1.2.landmarks[1:24,,]
Set1.2.landmarks
dim(Set1.2.landmarks)

# Set 3 Rep 2
Scale.3.rep2<-Rep3.18[c(1:18)*52]
Scale.3.rep2
scale.3.rep2<-sub("SCALE=", "", Scale.3.rep2)
scale.3.rep2
mode(scale.3.rep2)<-"numeric"
is.numeric(scale.3.rep2)

# generate an empty array to store the 18 configurations (specimens) for Rep 2 of Test II
n<-18
Set3.2.land<-array(NA, dim=c(24,2,18)) 
for (i in 1:n){Set3.2.land[,,i]<-matrix(Photo2.18[(i*52-50):(i*52-3)],24,2,byrow=T)}

Set3.2.land[,,12]
mode(Set3.2.land)<-"numeric"

# set the scale
Set3.2.landmarks<-Set3.2.land
for (i in 1:18)
{Set3.2.landmarks[,,i]<-(Set3.2.land[,,i]*(scale.3.rep2[i]))}
Set3.2.landmarks<-Set3.2.landmarks[1:24,,]
Set3.2.landmarks
dim(Set3.2.landmarks)

# Same Specimen, Two different images
p3<-36
Reps2photo<-array(NA, dim=c(24,2,p3))
Reps2photo[,,1:18]<-Set1.2.landmarks
Reps2photo[,,19:36]<-Set3.2.landmarks
dim(Reps2photo)
plot(Reps2photo)

Rep2sup<-gpagen(Reps2photo)

indphoto1<-1:18
indphoto2<-1:18
indsPhotos<-c(indphoto1, indphoto2)
indsPhotos<-as.factor(indsPhotos)

# model 3, size by individuals
size2<-Rep2sup$Csize
mod3<-summary(aov(size2~indsPhotos))
s2within_rep2<- mod3[[1]][2,3]
s2among_rep2<- (mod3[[1]][1,3] - s2within)/2
#We compute the % of measurement error for size.
s2within_rep2/(s2within_rep2+s2among_rep2) # 3% error in size by photo

shape2<-Rep2sup$coords
S2<-two.d.array(shape2)
#For shape, things are similar, we must estimate SS and crossproduct matrices first.
# model 4 shape by individuals in different photos
mod4<-lm(S2~indsPhotos)
SSres2<-crossprod(residuals(mod4)) 
SSeffect2<-crossprod(scale(mod4$fitted, scale=F))
tracewithin2<-sum(diag(SSres2))/ mod4$df.residual
traceamong2<-(sum(diag(SSeffect2))/(length(unique(indsPhotos))-1) - tracewithin2)/2
#We compute the % of shape measurement error.
tracewithin2/(tracewithin2+traceamong2) # 11.9% shape measurement error by image

