### Supplementary Material of the manuscript: "Porosity, openness, and exposure: 
### Identification of underlying factors associated with semi-outdoor spaces’ thermal performance and clustering in tropical high-density Singapore".
# written by Juan Gamero-Salinas, June 2022



# Loading dataframe
setwd("G:/...")
sos <- read.table("sos.txt", header=TRUE)
attach(sos)


# Preparing data for EFA
# Load libraries for EFA
library(factoextra)
library(corrplot)
library(psych)
library(GPArotation)


# Subset dataframe with only spatial attributes
df <- sos[,8:20] # 13 spatial attributes
df


### Exploratory Factor Analysis (EFA)

## Performing preliminary tests

# Pearson's Correlation Matrix
head(df,63)
df_matrix <- corr.test(df)
df_matrix 
str(df_matrix)
write.csv(df_matrix$ci,"correlationattributes.csv") # Save matrix in .csv file
corrplot(cor(df),method = 'circle', type = c("lower"), addCoef.col = 1,  number.cex = 0.6, tl.cex = 0.75, tl.col = 1, number.font = 1) # Correlogram
corrplot(cor(df),method = 'circle', type = c("lower"), addCoef.col = 1,  number.cex = 0.6, tl.cex = 0.75, tl.col = 1, number.font = 1, order ="hclust")


# Verifying the degree of adequacy of the factor analysis: KMO's MSA
KMO(df) 

# Verifying the degree of adequacy of the factor analysis: Bartlett's Test of Sphericity
cor.matrix <- cor(df)
cortest.bartlett(cor.matrix, n = nrow(df))



## Determining the estimation method for factor analysis

# Checking kurtosis and skewness of variables
kurtosi(df)
skew(df)

# Log-transforming variables outside range of -1.5 & +1.5
df$perimeter = log(df$perimeter)
df$totalfrontage = log(df$totalfrontage)
df$openfrontage = log((df$openfrontage + 1)) # To avoid NaNs due to negative values
df$volume = log(df$volume)
df$area = log(df$area)
df$solid = log(df$solid)

# Checking again kurtosis and skewness of variables
kurtosi(df)
skew(df)





## Determining the number of factors 

# Parallel Analysis
parallel <- fa.parallel(df, fm="ml", fa="fa", show.legend = FALSE)
parallel
str(fa.parallel)
str(parallel)




## Rotating and interpreting the factor pattern matrix

# Selecting rotation method: Oblique rotation
fa.mydata <- fa(r=df, nfactors = 3, rotate = "promax", fm ="ml", scores = "Bartlett")
print(fa.mydata, cut = 0.40) # Attributes

# Orientation and HFG need to be removed from the analysis since they have factor loadings < 0.4.
# Parallel Analysis & Rotation has to be performed again without those variables

df$Orientation <- NULL    # Removing Orientation 
df$HFG <- NULL            # Removing HFG
head(df,63)

# Running Parallel Analysis again
parallel <- fa.parallel(df, fm="ml", fa="fa", show.legend = FALSE) # Result: 3 factors
parallel # Print eigenvalues
str(fa.parallel)
str(parallel)

# Executing Oblique Rotation again
fa.mydata <- fa(r=df, nfactors = 3, rotate = "promax", fm ="ml", scores = "Bartlett")
print(fa.mydata, cut = 0.40) # Attributes

# Results:
# ML1 (Volume porosity): void + height + solid + volume
# ML2 (Perimeter openness): totalfrontage + openfrontage + perimeter
# ML3 (Exposure to sky): svf + gnpr + area + depth

# Estimation of factor scores
head(fa.mydata$scores, 63)
round(cor(fa.mydata$scores,use="pairwise"),3)

# Adding factor scores columns to original dataframe
scores = fa.mydata$scores
scores


  
  
### Multivariate regression analysis
  

# Convert to dataframe
as.data.frame(scores)
df_scores <- as.data.frame(scores)
df_scores



# Append scores to dataframe
sos$vporosity <- df_scores$ML1
sos$popeness <- df_scores$ML2
sos$exposure <- df_scores$ML3
head(sos)

# Load packages to calculate standardise coefficients and Condition Index (CI)
library(misty)
library(lm.beta)


## Perform multivariate regressions: validating retained factors
model.ta <- lm(Ta ~ vporosity + popeness + exposure, data=sos) # Ta
summary(model.ta)
coef_lmbeta <- lm.beta(model.ta)
coef_lmbeta
collin.diag(model.ta, print = c("all"), digits = 3, p.digits = 3,
            check = TRUE, output = TRUE)

model.tmrt <- lm(Tmrt ~ vporosity + popeness + exposure, data=sos) # Tmrt
summary(model.tmrt)
coef_lmbeta <- lm.beta(model.tmrt)
coef_lmbeta
collin.diag(model.tmrt, print = c("all"), digits = 3, p.digits = 3,
            check = TRUE, output = TRUE)

model.va <- lm(Va ~ vporosity + popeness + exposure, data=sos) # Va
summary(model.va)
coef_lmbeta <- lm.beta(model.va)
coef_lmbeta
collin.diag(model.va , print = c("all"), digits = 3, p.digits = 3,
            check = TRUE, output = TRUE)

model.rh <- lm(RH ~ vporosity + popeness + exposure, data=sos) # RH
summary(model.rh)
coef_lmbeta <- lm.beta(model.rh)
coef_lmbeta
collin.diag(model.rh , print = c("all"), digits = 3, p.digits = 3,
            check = TRUE, output = TRUE)

model.pmv1 <- lm(pmvg1 ~ vporosity + popeness + exposure, data=sos) # Gagge's PMV for 1 MET
summary(model.pmv1)
coef_lmbeta <- lm.beta(model.pmv1)
coef_lmbeta
collin.diag(model.pmv1, print = c("all"), digits = 3, p.digits = 3,
            check = TRUE, output = TRUE)

model.pmv1.5 <- lm(pmvg1.5 ~ vporosity + popeness + exposure, data=sos) # Gagge's PMV for 1.5 METs
summary(model.pmv1.5)
coef_lmbeta <- lm.beta(model.pmv1.5)
coef_lmbeta
collin.diag(model.pmv1.5, print = c("all"), digits = 3, p.digits = 3, 
            check = TRUE, output = TRUE)

model.pmv2 <- lm(pmvg2 ~ vporosity + popeness + exposure, data=sos) # Gagge's PMV for 2 METs
summary(model.pmv2)
coef_lmbeta <- lm.beta(model.pmv2)
coef_lmbeta
collin.diag(model.pmv2 , print = c("all"), digits = 3, p.digits = 3,
            check = TRUE, output = TRUE)

model.set1 <- lm(set1 ~ vporosity + popeness + exposure, data=sos) # Gagge's SET for 1 MET
summary(model.set1)
coef_lmbeta <- lm.beta(model.set1)
coef_lmbeta
collin.diag(model.set1, print = c("all"), digits = 3, p.digits = 3,
            check = TRUE, output = TRUE)

model.set1.5 <- lm(set1.5 ~ vporosity + popeness + exposure, data=sos) # Gagge's SET for 1.5 METs
summary(model.set1.5)
coef_lmbeta <- lm.beta(model.set1.5)
coef_lmbeta
collin.diag(model.set1.5, print = c("all"), digits = 3, p.digits = 3, 
            check = TRUE, output = TRUE)

model.set2 <- lm(set2 ~ vporosity + popeness + exposure, data=sos) # Gagge's SET for 2 METs
summary(model.set2)
coef_lmbeta <- lm.beta(model.set2)
coef_lmbeta
collin.diag(model.set2 , print = c("all"), digits = 3, p.digits = 3,
            check = TRUE, output = TRUE)

model.pet1 <- lm(pet1 ~ vporosity + popeness + exposure, data=sos) # PET for 1 MET
summary(model.pet1)
coef_lmbeta <- lm.beta(model.pet1)
coef_lmbeta
collin.diag(model.pet1, print = c("all"), digits = 3, p.digits = 3,
            check = TRUE, output = TRUE)

model.pet1.5 <- lm(pet1.5 ~ vporosity + popeness + exposure, data=sos) # PET for 1.5 METs
summary(model.pet1.5)
coef_lmbeta <- lm.beta(model.pet1.5)
coef_lmbeta
collin.diag(model.pet1.5, print = c("all"), digits = 3, p.digits = 3, 
            check = TRUE, output = TRUE)

model.pet2 <- lm(pet2 ~ vporosity + popeness + exposure, data=sos) # PET for 2 METs
summary(model.pet2)
coef_lmbeta <- lm.beta(model.pet2)
coef_lmbeta
collin.diag(model.pet2 , print = c("all"), digits = 3, p.digits = 3,
            check = TRUE, output = TRUE)



  
### Hierarchical Clustering with agnes() function in cluster package

# Select only factors with their factor scores
sos_clustering <- sos[,c("vporosity","popeness", "exposure")]
sos_clustering

# Load library for performing hierarchical clustering
library(cluster) # Load cluster() package


## Validating factor scores for clustering
# Define linkage method
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# Function to compute agglomerative coefficient
ac <- function(x) {agnes(sos_clustering, method = x)$ac} # calculate agglomerative coefficient

# Calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac)

# We can see that Ward's minimum variance method produces the highest 
# agglomerative coefficient

# Perform hierarchical clustering using Ward's minimum variance
clust <- agnes(sos_clustering, method = "ward")
pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram") 





# Determine the Optimal Number of Clusters

# Calculate the number of clusters
library(NbClust)
NbClust(data = sos_clustering, diss = NULL, distance = "euclidean",
        min.nc = 2, max.nc = 8, method = "ward.D2")
# According to majority rule: 3 clusters


#compute distance matrix
library(dplyr)
d <- dist(sos_clustering, method = 'euclidean')
final_clust <- hclust(d, method = 'ward.D2') # ward in agnes() function corresponds to ward.D2 in hclust() function
final_clust

# Plot and customize dendogram
dend <- as.dendrogram(final_clust)
plot(dend)

library(dendextend)
dend %>% set("branches_k_color", value = c("#FF0000","#002060","#70AD47"), k = 3) %>% set("labels_cex", 0.5) %>% set("branches_lwd", 2) %>% plot(horiz = TRUE)

# Add clusters to original dataframe
groups <- cutree(final_clust, k = 3)
table(groups)
groups

rect.hclust(final_clust , k = 3, border = 2:6)
sos_factors_clusters <- cbind(sos, cluster = groups)
sos_factors_clusters # prints dataframe with factor scores and clusters





# Generate mean values of retained factors for each cluster
aggregate(sos_factors_clusters, by=list(cluster=sos_factors_clusters$cluster), mean)
write.xlsx(head(sos_factors_clusters, 63),"sos_factors_clusters.xlsx") # Export it to .csv





### Visualization of factor scores for each cluster

# Load library
library(plotly)

# Assign columns to x, y & z axis
x <- sos_factors_clusters$vporosity
y <- sos_factors_clusters$popeness
z <- sos_factors_clusters$exposure

# Create legend
sos_factors_clusters$cluster[which(sos_factors_clusters$cluster == 1)] <- 'Horizontal Breezeway (HB)'
sos_factors_clusters$cluster[which(sos_factors_clusters$cluster == 2)] <- 'Perimeter Buffer (PB)'
sos_factors_clusters$cluster[which(sos_factors_clusters$cluster == 3)] <- 'Vertical Breezeway (VB)'
sos_factors_clusters$cluster <- as.factor(sos_factors_clusters$cluster)

# Execute code
fig <- plot_ly(sos_factors_clusters, x = ~x, y = ~y, z = ~z, color = ~cluster, colors = c("#70AD47","#FF0000","#002060"), opacity = 0.5)
fig <- fig %>% add_markers()          
fig <- fig %>% layout(scene = list(xaxis = list(title = 'VP'),
                                   yaxis = list(title = 'PO'),
                                   zaxis = list(title = 'ES')))

fig # print




