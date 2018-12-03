#Project-I
#Data Mining in Petroleum Engineering
#Pritesh Bhoumick
#OU ID - 113292019
#-------------------------------------

rm(list=ls(all=TRUE)) #Clear all existing variables and data entries
install.packages("xlsx")
library(xlsx)
install.packages("kohonen") #Install Kohonen SOM package
library(kohonen) #Attach Kohonen SOM method library

# ----------------------- READ DATA FROM INPUT FILE --------------------------

Alldata=read.delim("Input_data.las", sep="\t", header = TRUE) #Read input data from default directory
attach(Alldata) #Attach data from input file

# ------------------------ SOM MAPPING WITH 10x10 GRID --------------------------

Inputdata_SOM=Alldata[,3:7] #Extract relevant columns for analysis

#Convert data frame to matrix and center and scale variables
Inputdata_scale_SOM=as.matrix(scale(Inputdata_SOM)) 

set.seed(5) #Set seed value to direct algorithm to follow constant path

#Create SOM Grid and train grid prior to SOM
SOM_grid=somgrid(xdim = 10, ydim = 10, topo = "hexagonal") 

#Train SOM for no. of iterations, learing rates and neighbourhood available 
SOM_model=som(Inputdata_scale_SOM,grid=SOM_grid,rlen=300,
              alpha=c(0.05,0.01),keep.data=TRUE,n.hood="circular") 

#Plot iteration training progress with time
plot(SOM_model,type = "changes")
#Plot no. of node counts map
plot(SOM_model,type = "count") 
#Plot visualization of the distance between each node and its neighbours
plot(SOM_model,type = "dist.neighbours")
#Plot normalized values of original variables to generate SOM
plot(SOM_model,type = "codes") 

par(mfrow=c(3,3)) #Create plot map structure
plot(SOM_model, type = "property", property = SOM_model$codes[,1], 
     main = names(Inputdata_SOM)[1]) #Property map for RHOB
plot(SOM_model, type = "property", property = SOM_model$codes[,2], 
     main = names(Inputdata_SOM)[2]) #Property map for PEF
plot(SOM_model, type = "property", property = SOM_model$codes[,3], 
     main = names(Inputdata_SOM)[3]) #Property map for NPHI
plot(SOM_model, type = "property", property = SOM_model$codes[,4], 
     main = names(Inputdata_SOM)[4]) #Property map for GR
plot(SOM_model, type = "property", property = SOM_model$codes[,5], 
     main = names(Inputdata_SOM)[5]) #Property map for AT90

#Using hierarchical clustering to cluster the data
SOM_cluster = cutree(hclust(dist(SOM_model$codes)), 3)

#Define Palette color coding for clustered map
color_palette = c('#d62728', '#a306d2', '#4356d2','#e376c2')

#Show clustered plots with color coding
plot(SOM_model, type="codes", bgcol = color_palette[SOM_cluster], 
     main = "Clusters")
add.cluster.boundaries(SOM_model, SOM_cluster)

#Finding out the optimal no. of facies from the input data using SSB and SSW
data_SOM = SOM_model$codes
ssw = (nrow(data_SOM)-1)*sum(apply(data_SOM,2,var))
ssb = (nrow(data_SOM)-1)*sum(apply(data_SOM,2,var))
for (i in 1:15) ssw[i] = sum(kmeans(data_SOM, centers=i)$withinss)
for (i in 1:15) ssb[i] = sum(kmeans(data_SOM, centers=i)$betweenss)

#Plotting SSW and SSB
plot(1:15, ssb, type="b", col="red", xlab="Number of Clusters", ylab = "",
     main="Cluster sum of squares", ylim = range(ssb,ssw))
lines(1:15, ssw, type="b", col="blue")
legend("right",col=c("red","blue"),lty=1,legend=c("SSB","SSW"))

# ----------------------- SOM MAPPING WITH 3x1 MAPPING -----------------------

Inputdata_SOM2=Alldata[,3:7] #Extract relevant columns for analysis

#Convert data frame to matrix and center and scale variables
Inputdata_scale_SOM2=as.matrix(scale(Inputdata_SOM2)) 

#Set seed value to direct algorithm to follow constant path
set.seed(5)

#Create SOM Grid and train grid prior to SOM
SOM_grid_=somgrid(xdim = 3, ydim = 1, topo = "hexagonal") 

#Train SOM for no. of iterations, learing rates and neighbourhood available 
SOM_model_=som(Inputdata_scale_SOM2,grid=SOM_grid_,rlen=300,
              alpha=c(0.05,0.01),keep.data=TRUE,n.hood="circular") 

#Using hierarchical clustering to cluster the data
SOM_cluster_ = cutree(hclust(dist(SOM_model_$codes)), 3)

#Define Palette color coding for clustered map
color_palette = c('#d62728', '#a306d2', '#4356d2','#e376c2')

#Show clustered plots with color coding
plot(SOM_model_, type="codes", bgcol = color_palette[SOM_cluster_], 
     main = "Clusters")
add.cluster.boundaries(SOM_model_, SOM_cluster_)

#Extract the facies from the clustered data and save it in a excel file as datasheet
Unit = map(SOM_model_, Inputdata_scale_SOM2)
Facies_SOM = as.matrix(Unit$unit.classif)
WELL = Alldata$WELL
DEPTH = Alldata$DEPTH
FACIES = Facies_SOM[,1]
Facies_SOM = cbind.data.frame(WELL, DEPTH, FACIES)
colnames(Facies_SOM) = c("Well", "Depth", "Facies") #Assign column name to the data
write.xlsx(Facies_SOM, "Facies_SOM.xlsx")

# ----------------------------- K-MEANS CLUSTERING ----------------------------------

Inputdata_km=Alldata[,3:7] #Extract relevant columns for analysis

#Convert data frame to matrix and center and scale variables
Inputdata_scale_km=as.matrix(scale(Inputdata_km)) 

#Set seed value to direct algorithm to follow constant path
set.seed(5)

#Finding out the optimal no. of facies from the input data using SSB and SSW
kmeans_result = kmeans(Inputdata_scale_km,100)
data_km = kmeans_result$centers
ssw = (nrow(data_km)-1)*sum(apply(data_km,2,var))
ssb = (nrow(data_km)-1)*sum(apply(data_km,2,var))
for (i in 1:15) ssw[i] = sum(kmeans(data_km, centers=i)$withinss)
for (i in 1:15) ssb[i] = sum(kmeans(data_km, centers=i)$betweenss)

#Plotting SSW and SSB
plot(1:15, ssb, type="b", col="red", xlab="Number of Clusters", ylab = "",
     main="Cluster sum of squares", ylim = range(ssb,ssw))
lines(1:15, ssw, type="b", col="blue")
legend("right",col=c("red","blue"),lty=1,legend=c("SSB","SSW"))

#Apply Principal Component Analysis method
pca = princomp(Inputdata_scale_km, cor = T)
pca_comp = pca$scores

#Calculate two principal components in PCA
PCA_Component1 = pca_comp[,1]
PCA_Component2 = pca_comp[,2]
kmeans_analysis = cbind(PCA_Component1,PCA_Component2)

#Create clusters based on K-Mans algorithm
kcluster = kmeans(kmeans_analysis,3)

#Plot the two principal components to see the clusters
plot(PCA_Component1,PCA_Component2,col=kcluster$cluster)
points(kcluster$centers,pch=9, col="white", cex=1)

#Extract the facies from the clustered data and save it in a excel file as datasheet
Facies_km = as.matrix(kcluster$cluster)
WELL = Alldata$WELL
DEPTH = Alldata$DEPTH
FACIES = Facies_km[,1]
Facies_km = cbind.data.frame(WELL, DEPTH, FACIES)
colnames(Facies_km) = c("Well", "Depth", "Facies") #Assign column name to the data
write.xlsx(Facies_km, "Facies_Kmeans.xlsx")