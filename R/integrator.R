## Description

#Our standard network contains 365 drugs in it, with 8 of those drugs having to be
#imputed in the perturbation layer since we're using the new perturbation signatures.
#The aim of this notebook is to create the largest network possible by using all information
#from the sensitivity layer, and not subsetting to just the drugs that are in common with
#structure and perturbation.

# source('~/Desktop/DNFAPI/R/structureDataFlexible.R')
# source('~/Desktop/DNFAPI/R/snfModifiedHelpers.R')
# source('~/Desktop/DNFAPI/R/sensitivityDataFlexible.R')
# source('~/Desktop/DNFAPI/R/perturbationDataFlexibleCustomSig.R')
# source('~/Desktop/DNFAPI/R/medianSimilarity.R')

integrator <- function(struct, perturb, sensi){
  ## Full Network on Combined Benchmark

  #Get all necessary imports and functions
  set.seed(9833)

  badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

  #Load LINCS metadata file. Used it contains the smiles string for the drugs in the perturbation layer. These are then used to construct the structure layer.
  lincs.meta <- rbind(struct80, struct)
  lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
  lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

  #Load in senstivity, imaging data as well as perturbation signatures

  sens.data <- sens
  pert.data <- cbind(pert80, perturb) #pert80 # #PerturbationDataFlexibleCustomSig(pert.file.name)  ## 978 X 23

  #Subset the lincs metadata file to only the drugs appearing in the sensitivity layer

  sens.names <- rownames(sens.data)
  #return(list(sens=sens.names, struc=lincs.meta$pert_iname))
  lincs.meta.subset <- lincs.meta[match(sens.names, lincs.meta$pert_iname),]
  lincs.meta.subset <- lincs.meta.subset[!is.na(lincs.meta.subset$X),]

  #Now we create the structure fingerprints based on the subsetted lincs metadata

  strc.data <- strc.data #StructureDataFlexible(lincs.meta.subset)  ## a vector  --> 239 elemnts

  ## Find Common Drugs

  #Now to find the drugs that are common between all layers

  #Now we can choose the layers to be used in the integration and see the drugs that are common between them. We leave out the drug names from the perturbation layer so common.drugs will be a vector of 365 drug names.

  layers <- list(sens.names = sort(colnames(sens.data)),
                 strc.names = names(strc.data))
  common.drugs <- Reduce(intersect, Filter(Negate(is.null),layers))

  ## Create Similarities

  #At this point we have all the data necessary to create the various similarity matrices.
  #Keep in mind that in the case of the senstivity layer, the data loaded in is already a similarity matrix, so there is nothing to do there.
  strc.cor <- fingerprint::fp.sim.matrix(strc.data, method = "tanimoto")
  colnames(strc.cor) <- names(strc.data)
  rownames(strc.cor) <- names(strc.data)

  strc.cor <- strc.cor[common.drugs, common.drugs]

  pert.cor <- cor(pert.data[, intersect(common.drugs, colnames(pert.data))], method="pearson", use="pairwise.complete.obs")

  pert.cor <- medianSimilarity(list(pert.cor))[[1]]

  sens.cor <- sens.data[common.drugs, common.drugs]

  ## Integrate Network

  #all.drugs is the list of 365 drug names. CreateAugmentedMatrixSkeletons() will create a brand new matrix for each layer, and since the pertubation layer is missing 8 drugs relative to sensitivity and structure, the new skeleton matrix for this layer will contain these 8 drugs, whose values will be imputed with values from the sensitivity and structure layers.

  all.drugs <- sort(rownames(sens.cor))

  correlation.matrices <- list(sens=sens.cor, pert=pert.cor, strc=strc.cor)

  affinity.matrices <- CreateAffinityMatrices(correlation.matrices)
  augmented.matrices <- CreateAugmentedMatrixSkeletons(names(correlation.matrices), all.drugs)
  augmented.matrices <- ReplaceAugmentedExistingValues(augmented.matrices, affinity.matrices)
  affinity.matrices <- ReplaceAffinityMatrixValuesFast(augmented.matrices, correlation.matrices,
                                                       all.drugs)
  affinity.matrices <- medianSimilarity(affinity.matrices)

  integrated <- SNFtool::SNF(affinity.matrices)
  rownames(integrated) <- all.drugs
  colnames(integrated) <- all.drugs
  return(integrated)
}
