# Easy way to convert sociomatrix to edgelist

library(statnet)

filenames_sm <- list.files(path = "./normalized_sociomatrices")

for (x in filenames_sm) {
	t <- read.csv(paste("./normalized_sociomatrices/", x, sep = ""), header = FALSE)
	e <- as.edgelist.sna(data.matrix(t))
	write.table(e, paste("./edgelists/e_", x, sep = ""), sep = ",", row.names = FALSE, col.names = FALSE)
}


# Quantify network metrics

library(tnet)
library(igraph)
library(sna)
library(psych)

filenames_e <- list.files(path = "./edgelists")

write(c("filename", "n", "Newman modularity", "weighted distance", "weighted diameter", "weighted clustering", "eigenvector centralization"), "fulldata.txt", ncolumns = 7, append = FALSE, sep = "\t")

for (y in filenames_e) {
	print(y)

	a <- read.csv(paste("./edgelists/", y, sep = ""), header = FALSE)
	i <- a$V1; j <- a$V2; w <- a$V3; b <- cbind(i, j, w)

	ig <- graph.data.frame(b, directed = FALSE)
	dens <- graph.density(ig) / 2
	cent <- igraph::evcent(ig, weights = w, scale = FALSE)
	Ce.ei <- (-sum(cent$vector - max(cent$vector))) / sum(rep(1, cent$options$n) - cent$vector)
	Newman.mod <- leading.eigenvector.community(ig, weights = w)$modularity
	diam <- diameter(ig, directed = FALSE, weights = w)
	size <- cent$options$n

	c <- symmetrise_w(b)
	Cl.wgm <- clustering_w(c, measure = "gm")
	dist.mat <- distance_w(c)
	dist.w <- mean(dist.mat, na.rm = TRUE)
	# harmonic.mean(dist.mat)

	write(c(y, size, Newman.mod, dist.w, diam, Cl.wgm, Ce.ei), "fulldata.txt", ncolumns = 7, append = TRUE, sep = "\t")
}


# Run PGLS models on all results

library(geiger)
library(caper)
library(plyr)

iterations = 1100
twomulti_data <- read.csv("reduceddata.csv")
simpletree<-read.nexus("binindaemonds-revised-pruned.nex")


################
# centralization

# change filename to 'output_....txt'
write(c("iteration","beta","std error","t-value","p-value","R-squared","adjusted R-squared","lambda"), "output-centralization.txt", ncolumns=8, append=FALSE, sep="\t")

for (i in 1:iterations)
{
	subd_data <- ddply(twomulti_data, .(species), function(x) {x[sample(nrow(x), 1),]})
	
	row.names(subd_data) <- subd_data$species

	primate <- comparative.data(phy = simpletree, data = subd_data, names.col = species, vcv = TRUE)
	# change 'cent' to 'clust' 'mod' 'dist'
	try({model.pgls<-pgls(n ~ cent, data = primate, lambda='ML')
	
	# change filename to 'output_....txt'
	write(c(i,summary(model.pgls)$coefficients[2,1],summary(model.pgls)$coefficients[2,2],summary(model.pgls)$coefficients[2,3],summary(model.pgls)$coefficients[2,4],summary(model.pgls)$r.squared,summary(model.pgls)$adj.r.squared,summary(model.pgls)$param[2]), "output-centralization.txt", ncolumns=8, append=TRUE, sep="\t")})
}
