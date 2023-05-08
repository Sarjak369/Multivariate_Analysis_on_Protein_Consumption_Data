# MVA Midterm Exam
# Name: Sarjak Maniar
# Email: sm2732@scarletmail.rutgers.edu

library(factoextra)
library(corrplot)
library(cluster)
library(psych)
library(dplyr)
library(devtools)
library(readr)
library(psych)


# The accompanied dataset gives estimates of the average protein consumption in grams per 
# person per day) from different food sources for the inhabitants of 25 European countries.

## Importing Data and initial analysis
## Importing csv file

Protein_Consumption <- read.csv("/Users/sarju/Desktop/MITA Sem 2/MVA/MVA_Midterm/Protein_Consumption.csv")
Protein_Consumption
Protein_Consumption <- as.data.frame(Protein_Consumption)
Protein_Consumption

#Dimension of the dataset
dim(Protein_Consumption) # 25 11 [There are 25 rows and 11 columns]
colnames(Protein_Consumption) # printing all the column names

attach(Protein_Consumption)

# Printing first 5 rows of the dataset
head(Protein_Consumption)

# Printing the summary of the dataset
summary(Protein_Consumption)

# Printing the structure of the dataset
str(Protein_Consumption)
# As we can see that all the columns are of type int except the column Country, which is of type chr

# =============================================================================================

# Question 1 - Use principal components analysis to investigate the relationships between the 
# countries on the basis of these variables

# Getting the Correlations 
cor(Protein_Consumption[-1])

# Plotting correlation plot to understand the most correlated variables
Protein_Consumption_plot <- cor(Protein_Consumption[-1])
corrplot(Protein_Consumption_plot, method="circle")

# Inferences: 
# We can see from the correlation plot, that Cereal and Eggs are negatively correlated. 
# Also, White Meat and Eggs are positively correlated.

# Finding the principal components of data

# Using prcomp to compute the principal components (eigenvalues and eigenvectors). 
# With scale=TRUE, variable means are set to zero, and variances set to one
Protein_Consumption_pca <- prcomp(Protein_Consumption[,-1],scale=TRUE) 
#Scaling to Standardize the data values
Protein_Consumption_pca

# Sample scores -> stored in Protein_Consumption_pca$x 
Protein_Consumption_pca$x  # returns the sample scores for each country along each principal component.

# Sample scores (also called factor scores) are the values of each observation 
# (in this case, each country) along each principal component (PC).
# Sample scores basically tells how much each observation contributes to each PC

# Singular values (square roots of eigenvalues) -> stored in Protein_Consumption_pca$sdev
Protein_Consumption_pca$sdev

# Loadings (eigenvectors) -> stored in Protein_Consumption_pca$rotation
Protein_Consumption_pca$rotation # provides the matrix of loadings or principal component coefficients.

# Variable means -> stored in Protein_Consumption_pca$center
Protein_Consumption_pca$center

# Variable standard deviations -> stored in Protein_Consumption_pca$scale
Protein_Consumption_pca$scale

# A table containing eigenvalues and %'s accounted, follows
# Eigenvalues are sdev^2

#Extract variance against features
eigenvalues<-Protein_Consumption_pca$sdev^2
eigenvalues
sum(eigenvalues)
names(eigenvalues) <- paste("PC",1:10,sep="")
eigenvalues
sumoflambdas <- sum(eigenvalues)
sumoflambdas


#Variance %
bfcvar<- (eigenvalues/sumoflambdas)*100
round(bfcvar,10)
barplot(bfcvar,main="Scree plot",xlab="Principal Component",ylab="Percent Variation")

#Calculate cumulative of variance
cumvar <- cumsum(bfcvar)
cumvar
matlambdas <- rbind(eigenvalues,bfcvar,cumvar)
round(matlambdas,10)

eigenvec_Protein_Consumption <- Protein_Consumption_pca$rotation
eigenvec_Protein_Consumption

# Relation between the variables w.r.t countries
Protein_Consumption_env <- cbind(data.frame(Country),Protein_Consumption_pca$x)
Protein_Consumption_env

#Visualize PCA using Scree plot
screeplot(Protein_Consumption_pca, type='bar',main='Scree plot')
summary(Protein_Consumption_pca)


# Principal Component Analysis is used to reduce the large set of variables to a 
# small set that still contains most of the information from the large dataset and thereby, 
# reducing the dimension of the dataset.
# After performing PCA on the Protein Consumption Dataset, we obtain 10 Principal components. 
# Each of these components represent the percentage of variability present in the dataset.
# In other words, PC1 explains 41.3% of total variance, PC2 explains 17.4% and so on. 
# We will consider the first six principal components as they sum up to 93% of total variance 
# and the other four can be discarded as they contribute to only 7% of total variance.


# Better Ways to Visualize

library(factoextra)
library(FactoMineR)
library(ggfortify)
library(psych)
library(corrplot)
library(devtools)

# Correlation
pairs.panels(Protein_Consumption[,-1],
             gap = 0,
             bg = c("red", "blue")[Protein_Consumption$Country],
             pch=21)

pairs.panels(Protein_Consumption_pca$x,
             gap=0,
             bg = c("red", "blue")[Protein_Consumption$Country],
             pch=21)

# Variables - PCA (cos2)
fviz_eig(Protein_Consumption_pca, addlabels = TRUE)
fviz_pca_var(Protein_Consumption_pca,col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)

# The plot above is also known as Variable Correlation plots. 
# It shows the relationships between all variables. It can be interpreted as follow:
# Positively correlated variables are grouped together.
# Negatively correlated variables are positioned on opposite sides of the plot origin (opposed quadrants).
# The distance between variables and the origin measures the quality of the variables on the factor map. 
# Variables that are away from the origin are well represented on the factor map.


# Individuals - PCA (cos2)
fviz_pca_ind(Protein_Consumption_pca, col.ind = "cos2", 
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
             repel = TRUE)

# Biplot
biplot(Protein_Consumption_pca)

autoplot(Protein_Consumption_pca,
         data = Protein_Consumption[,-1],
         loadings = TRUE,
         labels = Protein_Consumption$Country)



# Different PCA Method. 
res.pca <- PCA(Protein_Consumption[,-1], graph = FALSE)
print(res.pca)

# Visualize and Interpret PCA using these functions 

# get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
# fviz_eig(res.pca): Visualize the eigenvalues
# get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
# fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
# fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.

eig.val <- get_eigenvalue(res.pca)
eig.val

# Scree Plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
# As we can see that we got the same plot with this method also.

var <- get_pca_var(res.pca)
# var$coord: coordinates of variables to create a scatter plot
# var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
# var$contrib: contains the contributions (in percentage) of the variables to the principal components. 
# The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).
var
# Coordinates
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)


# Correlation circle
fviz_pca_var(res.pca, col.var = "black")

# Quality of representation

corrplot(var$cos2, is.corr=FALSE)
# Total cos2 of variables on Dim.1 and Dim.2
# A high cos2 indicates a good representation of the variable on the principal component,
# -> In this case the variable is positioned close to the circumference of the correlation circle.

# A low cos2 indicates that the variable is not perfectly represented by the PCs,
# -> In this case the variable is close to the center of the circle.

fviz_cos2(res.pca, choice = "var", axes = 1:2)
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

# Change the transparency by cos2 values
fviz_pca_var(res.pca, alpha.var = "cos2")
corrplot(var$contrib, is.corr=FALSE)
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
fviz_pca_var(res.pca, alpha.var = "contrib")


# =============================================================================================

# Question 2 - Carry out cluster analysis to study relation between countries on their diet


# K-means Clustering

# Standardizing the data with scale()
matstd.Protein_Consumption <- scale(Protein_Consumption[-1])
matstd.Protein_Consumption

# kmeans() is the function for  non-hierarchical method. 
rownames(matstd.Protein_Consumption) <- Protein_Consumption$Country

# For 2 clusters, k-means = 2
kmeans2.Protein_Consumption <- kmeans(matstd.Protein_Consumption,2,nstart = 10)

# Computing the percentage of variation accounted for two  clusters
perc.var.2 <- round(100*(1 - kmeans2.Protein_Consumption$betweenss/kmeans2.Protein_Consumption$totss),1)
names(perc.var.2) <- "Perc. 2 clus"
perc.var.2

# For 3 clusters, k-means = 3
kmeans3.Protein_Consumption <- kmeans(matstd.Protein_Consumption,3,nstart = 10)
# Computing the percentage of variation accounted for three  clusters
perc.var.3 <- round(100*(1 - kmeans3.Protein_Consumption$betweenss/kmeans3.Protein_Consumption$totss),1)
names(perc.var.3) <- "Perc. 3 clus"
perc.var.3

# For 4 clusters, k-means = 4
kmeans4.Protein_Consumption <- kmeans(matstd.Protein_Consumption,4,nstart = 10)
# Computing the percentage of variation accounted for four clusters
perc.var.4 <- round(100*(1 - kmeans4.Protein_Consumption$betweenss/kmeans4.Protein_Consumption$totss),1)
names(perc.var.4) <- "Perc. 4 clus"
perc.var.4

# We divide the dataset into two clusters.
# Filtering properties which are in 1 cluster of k mean 2
clus1 <- matrix(names(kmeans2.Protein_Consumption$cluster[kmeans2.Protein_Consumption$cluster == 1]),
                ncol=1, nrow=length(kmeans2.Protein_Consumption$cluster[kmeans2.Protein_Consumption$cluster == 1]))
colnames(clus1) <- "Cluster 1"
clus1

# Filtering properties which are in 2 cluster of k mean 2
clus2 <- matrix(names(kmeans2.Protein_Consumption$cluster[kmeans2.Protein_Consumption$cluster == 2]),
                ncol=1, nrow=length(kmeans2.Protein_Consumption$cluster[kmeans2.Protein_Consumption$cluster == 2]))
colnames(clus2) <- "Cluster 2"
clus2




# We can use the hierarchical clustering method, which allows us to create a dendrogram 
# to visualize the relationships between the countries based on their protein consumption from different food sources.
# Hierarchical Clustering

# Calculating the distance matrix using Euclidean distance
dist.Protein_Consumption <- dist(Protein_Consumption, method = "euclidean") # Distance matrix

# Performing hierarchical clustering using complete linkage
fit <- hclust(dist.Protein_Consumption, method="complete")

#Plotting the dendogram

plot(fit, main="Dendrogram of European countries based on diet")

#Cutting the tree into 2 clusters
groups <- cutree(fit, k=2)

#Plotting dendogram with red borders around the 2 clusters
rect.hclust(fit, k=2, border="red")


# From the dendrogram, we can see that the countries can be grouped into two main clusters: 
# one containing mostly Northern and Western European countries, and another containing mostly 
# Southern and Eastern European countries.

# The countries in the second cluster (Northern and Western Europe) generally consume more meat, 
# fish, milk, and cereals, while the countries in the first cluster (Southern and Eastern Europe) 
# consume more fruits, vegetables, and pulses/nuts/oilseeds.

# Overall, the hierarchical clustering analysis confirms some of the dietary patterns and 
# differences that we observed in the earlier principal components analysis.


# Clustering is a method of grouping together a set of objects together in such a way 
# that the objects in one cluster is similar to the objects in the same cluster than 
# objects present in the different cluster.
# We form 2 clusters for the given dataset as seen in the dendogram. 
# This covers all the variance present in the dataset.


# ??????????????????????????


library(cluster)
library(factoextra)
library(magrittr)
library(NbClust)

# We will use agnes function as it allows us to select option for data standardization, the distance measure and clustering algorithm in one single function

(dist.Protein_Consumption <- agnes(Protein_Consumption, metric="euclidean", stand=TRUE, method = "single"))
#View(dist.Protein_Consumption)

# Description of cluster merging
dist.Protein_Consumption$merge

# Dendogram
plot(as.dendrogram(dist.Protein_Consumption), xlab= "Distance",xlim=c(8,0),
     horiz = TRUE,main="Dendrogram")

#Interactive Plots
#plot(agn.employ,ask=TRUE)
plot(dist.Protein_Consumption, which.plots=1)
plot(dist.Protein_Consumption, which.plots=2)
plot(dist.Protein_Consumption, which.plots=3)


# =============================================================================================

# Question 3 - Identify the important factors underlying the observed variables 
# and examine the relationships between the countries with respect to these factors


# Computing Correlation Matrix
corrm.Protein_Consumption <- cor(Protein_Consumption[c(-1,-11)])
corrm.Protein_Consumption
corrplot(corrm.Protein_Consumption,method="number")
fit.bfc <- principal(Protein_Consumption[c(-1,-11)], nfactors=4, rotate="varimax")
fit.bfc


# Note that here, we are excluding/ignoring the column 'total', because it is just having the summation of all the column values

# We are considering 4 factors, because of VSS (Very Simple Structure)
round(fit.bfc$values, 3)

# Loadings
fit.bfc$loadings
# Factor Loadings: This table displays the factor loadings for each variable. 
# It represents the correlation between each variable and each factor. 

# Communalities
fit.bfc$communality
# Communalities are estimates of the variance in each observed variable that can be 
# explained by the extracted factors.
# It shows the proportion of each variable's variance that can be explained by the 
# factors, e.g., Red Meat explains 0.862, White Meat explains 0.913, and so on.

# Rotated factor scores
head(fit.bfc$scores)
round(fit.bfc$values,3)

# Factor recommendation
fa.parallel(Protein_Consumption[c(-1,-11)])
# From this, we can inference that 4th component could be the best choice for number of factors

# Correlations within Factors
fa.plot(fit.bfc)

# Visualizing the relationship
fa.diagram(fit.bfc)

# Factor recommendations for a simple structure
vss(Protein_Consumption[c(-1,-11)]) 
summary(vss(Protein_Consumption[c(-1,-11)]) )

# Overall, the factor analysis suggests that there are four underlying factors that 
# explain most of the variance in the dataset. These factors represent different aspects 
# of the European diet, such as a diet high in Milk, Egg and Fish (factor 1), a diet high 
# in fish, starchy food and cereals (factor 2), a diet high in white meat, pulses/nuts/oilseeds  
# (factor 3), and a diet high in fruits and vegetables (factor 4). 
# The factor scores for each country provide a way to compare countries based on their diet patterns and 
# identify similarities and differences between them.


# Factor Analysis is the method of identifying the latent relational structure among 
# a set of variables and then narrowing it down to a smaller number of variables.
# We can see that we have reduced the factors to five, this contains most of the 
# information in the dataset.
# The factors RC1, RC2, RC3, RC4 can be helpful in analyzing the entire dataset.


















