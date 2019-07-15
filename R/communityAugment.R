communityAugment <- function(network, GMT_TARG) {

  ## all communities from integrative analysis results for nci60/l1000 or ctrpv2/l1000

  apcomb <- apcluster(network, q=0.9)
  ll <- list()
  for(i in 1:length(apcomb)){
    xx <- names(apcomb[[i]])
    ll[[i]] <- xx
  }

  indx <- sapply(ll, length)
  #indx <- lengths(lst)
  res <- as.data.frame(do.call(rbind,lapply(ll, `length<-`,max(indx))))
  llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
  row.names(llx) <- names(apcomb@exemplars) #the non-refined version
  result <- applyCommunity(network, llx)

  first <- llx

  dim(res)

  ll <- list()
  for(i in 1:length(GMT_TARG)){
    xx <- GMT_TARG[[i]]
    ll[[i]] <- xx
  }

  indx <- sapply(ll, length)
  #indx <- lengths(lst)
  res <- as.data.frame(do.call(rbind,lapply(ll, `length<-`,max(indx))))
  llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
  row.names(llx) <- unlist(names(GMT_TARG))
  #write.csv(llx," GMT_targ_chembl.csv", row.names=TRUE)
  dim(res)


  ## Keep communities with at least 2 drugs showing a known mechanism of action from GMT
  GMT_TARG2 <- c(as.character(unlist(GMT_TARG)))

  #Refine the Clusters: exclude drugs from each cluster which are not found in the GMT benchmark drugs
  Clust <- apcomb@clusters
  ClusterList <- lapply(Clust,function(x){if (length(intersect(GMT_TARG2 ,names(x[]))) >= 2) {names(x[])} else {NULL}})
  ClusterByCondition<- sapply(ClusterList, function(x) length(x) > 1)
  ClusterListRefined <- ClusterList[ClusterByCondition]

  indx <- sapply(ClusterListRefined, length)
  res <- as.data.frame(do.call(rbind,lapply(ClusterListRefined, `length<-`,max(indx))))
  llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res) #the refined version
  rownames(llx) <- rownames(first[match(llx$V1, first$V1), ])
  return(applyCommunity(result, llx))
}


#Add 1 to all drugs that interact with a representative community drug
applyCommunity <- function(network, communities){
  #finalDrugs <- c()
  for (row in 1:nrow(communities)) { #iterate through each row of the communities

    #Find the current community
    community <- communities[row, 3:ncol(communities)]
    community <- as.vector(droplevels(unlist(community)))

    #Find the representative of the community
    representative <- rownames(communities)[row]

    #Find the array positions of the drugs in the network that are also in this community
    drugPositions <- which(colnames(network) %in% community)

    #Augment all of these drug positions
    network[representative, drugPositions] <- network[representative, drugPositions] + 1

    # finalDrugs <- append(finalDrugs, drugPositions) #old line that was used for eliminating drugs that
  }
  #finalDrugs <- sort(unique(finalDrugs))
  #return(network[finalDrugs, finalDrugs])
  return(network)
}

