#Brian J. Spiesman and Claudio Gratton
#Flexible foraging null model analysis
#22 September 2015


#requires libraries: vegan, igraph, and bipartite

###Model input
#net = name of realized (local) weighted network, a matrix with plant names making up the row names and pollinator names making up the column names
#metanet = name of regional weighted network formed by combining all realized networks in the region, a matrix with plant names as row names and pollinator names as column names 
#phen = name of binary matrix indicating which plant and pollinator species have local temporal overlap in activity (same dimensions as net)
#n = number of randomizations to perform
#s.prob = must be either 1 (Null model 1: randomization weighted by marginal totals) or 2 (Null model 2: equal probability randomization).  

nullmod(net = net, metanet = metanet, phen = phen, n = 1000, s.prob = 1)

nullmod = function(net, metanet, phen, n, s.prob){
	require(vegan); require(igraph); require(bipartite)

	###Generate realized network values (all scaled 0-100)
	NODF = nestednodf(net, order=TRUE, weighted=T)$statistic[3]			#Realized Nestedness (NODF)

	g = graph.incidence(net, weighted=TRUE)
	wtc = max(walktrap.community(g)$modularity)
	fgc = max(fastgreedy.community(g)$modularity)
	Mod = 100*(max(wtc,fgc))								#Realized Modularity

	H2 = H2fun(round(net))[1]*100								#Realized Complementary Specialization (H2)

	bnet = decostand(net, method="pa")							#binary version of net
	Conn = 100 * sum(bnet)/(dim(bnet)[1]*dim(bnet)[2])				#Realized Connectance

	###Generate null models
	nullmods = vector("list", n)
	for (i in 1:n){
		pnet = metanet[rownames(net), colnames(net)] * phen 			#potential network: local species w/ regional interactions, corrected for local phenology mismatch
		bpnet = decostand(pnet, method="pa")					#binary verison of pnet
		p.r = rowSums(bpnet)/sum(bpnet)						#row sums of binary potential network
		p.c = colSums(bpnet)/sum(bpnet)						#colunm sums of binary potential network
		if (s.prob == 1){									#selection probability (Null model 1): randomization is weighted by marginal totals
			P = 1/(bpnet * sqrt(tcrossprod(p.r, p.c))); P[is.nan(P)] = 0; P[is.infinite(P)] = 0			
		}
		if (s.prob == 2){									#selection probability (Null model 2): randomization is not wieghted
			P = bpnet				
		}
		for (j in 1:ncol(net)){
			if (colSums(bpnet)[j] == 1){						#all species must have at least one partner - so if a pollinator has only 1 parnter its selection prob = 0
				P[,j] = 0 * bpnet[,j]
			}
		}
		for (j in 1:nrow(net)){
			if (rowSums(bpnet)[j] == 1){						#if a plant has only 1 parnter its selection prob = 0
				P[j,] = 0 * bpnet[j,]
			}
		}
		randomnet = bpnet									#set initial network for random deletion of interactions
		r_conn = 100 * sum(randomnet)/(nrow(randomnet)*ncol(randomnet))	#determine initial connectance of randomized network
		while (Conn < r_conn){								#loop through to pare down random binary network, stop when random and realized network connectance is equal
			sel = sample(1:length(net), 1, prob = P, replace = FALSE)
			Temp_randomnet = randomnet
			Temp_randomnet[sel] = 0
			check = min(rowSums(Temp_randomnet)) * min(colSums(Temp_randomnet))
				if (check >= 1){
					randomnet[sel] = 0
					P[sel] = 0
				}
			r_conn = 100 * sum(randomnet) / (nrow(randomnet)*ncol(randomnet))
		}
	
		p.r = rowSums(pnet)/sum(pnet)
		p.c = colSums(pnet)/sum(pnet)
		if (s.prob == 1){ 								#selection probability (Null model 1): randomization is weighted by marginal totals
			P = randomnet * sqrt(tcrossprod(p.r, p.c))				
		}
		if (s.prob == 2){									#selection probability (Null model 2): randomization is not wieghted
			P = bpnet				
		}
		n.int = round(sum(net)) 							#number of interactions in observed network
		n.int.randomnet = sum(randomnet)						#initial number of interactions in random network
		while (n.int.randomnet < n.int) {
			sel = sample(1:length(net), 1, prob = P, replace = TRUE)
			randomnet[sel] = 1 + randomnet[sel]		
			n.int.randomnet = sum(randomnet)
		}
		nullmods[[i]] = randomnet
	}

	###Null model analysis
	nullNODF = seq(1:n)

	#Nestedness
	for (i in 1:n){
		nNODF = nestednodf(nullmods[[i]], order=TRUE, weighted=TRUE)
		nullNODF[i] = nNODF$statistic[3]
	}
	p_NODF = length(nullNODF[nullNODF > NODF])/n					#Nestedness p-value
	z_NODF = (NODF -(mean(nullNODF)))/sd(nullNODF)					#Relative nestedness (z-value)

	#Modularity
	nullMod = seq(1:n)
	for (i in 1:n){
		nnet = as.matrix(nullmods[[i]])
		g = graph.incidence(nnet, weighted=TRUE)
		wtc = max(walktrap.community(g)$modularity)
		fgc = max(fastgreedy.community(g)$modularity)
		nullMod[i] = 100*(max(wtc,fgc))
	}
	p_Mod = length(nullMod[nullMod > Mod])/n						#Modularity p-value
	z_Mod = (Mod -(mean(nullMod)))/sd(nullMod)					#Relative modularity (z-value)

	#Specialization
	nullH2 = seq(1:n)										
	for (i in 1:n){
		nH2 = H2fun(nullmods[[i]])
		nullH2[i] = nH2[1]*100
	}
	p_H2 = length(nullH2[nullH2 > H2])/n						#Complementary specialization p-value
	z_H2 = (H2 -(mean(nullH2)))/sd(nullH2)						#Relative specialization (z-value)

	#Output for realized network values, p-values, and relative network values
	output = cbind(NODF, p_NODF, z_NODF, Mod, p_Mod, z_Mod, H2, p_H2, z_H2);row.names(output) = "model.output"
	print(output)
}
