map.paths = function(G.inferred, G){
	G.inferred = subGraph(intersect(nodes(G), nodes(G.inferred)), G.inferred)
	em = edgeMatrix(G.inferred)
	eval = data.frame()	
	for(i in 1:NCOL(em)){
		cand.path = paste(nodes(G.inferred)[em["from", i]], "->",nodes(G.inferred)[em["to", i]], sep="")
		path = sp.between(G, nodes(G.inferred)[em["from", i]], nodes(G.inferred)[em["to", i]])
		if(!is.na(path[[1]]$length)){
			eval = rbind(eval, data.frame(edge=cand.path, explained.by=paste(path[[1]]$path_detail, collapse="->")))
		}
		else
			eval = rbind(eval, data.frame(edge=cand.path, explained.by="none"))
	}
	eval
}

sp.path.compute = function(genelist, genelist2=NULL, G){
	require(RBGL)
	if(is.null(genelist2)){
		paths = c()		
		for(i in 1:length(genelist)){			
			for(j in 1:length(genelist)){
				cat(paste(i, "of", length(genelist), "*** and ***", j, "of", length(genelist), "\n",sep=" "))
				if(i != j){
					paths = c(paths, sp.between(G, genelist[i], genelist[j]))					
				}				
			}			
		}
		return(paths)				
	}	
	else{		
		paths = sp.between(G, genelist, genelist2)
		return(paths)
	}
}

plot.paths = function(paths, G, interactions, center=NULL, filename=NULL){
	require(igraph)
	if(!is.null(filename))
		pdf(filename)	
	nodes = unlist(strsplit(names(paths), ":"))
	nodes = union(nodes, unlist(sapply(paths, function(p) p$path_detail)))
	g = igraph.from.graphNEL(subGraph(setdiff(nodes, NA), G))
	vcol = rep("white", length(V(g)$name))
	if(!is.null(center))
		vcol[V(g)$name %in% center] = "red"
	e1 = get.edgelist(g)
	col = c()
	for(e in 1:NROW(e1)){
		take = which(as.character(interactions$Network.Object..FROM.) == e1[e,1] & as.character(interactions$Network.Object..TO.) == e1[e,2])
		if(length(take) > 0){
			if(interactions$Effect[take] == "Inhibition")
				col = c(col, "red")
			else if(interactions$Effect[take] == "Activation")
				col = c(col, "green")					
			else 
				col = c(col, "darkgrey")	
		}
	}
	tkplot(g, edge.color=col, vertex.color=vcol, vertex.label=V(g)$name, vertex.label.color="black", layout=layout.fruchterman.reingold, vertex.size=20)	
	if(!is.null(filename))
		dev.off()
	g
}

# example
#res = read.xlsx("MEF2CNetwork.xlsx", sheetIndex=2)
#G = new("graphNEL", edgemode="directed", nodes=union(as.character(net[,3]), as.character(net[,6])))
#G = addEdge(as.character(net[,3]), as.character(net[,6]), G)



