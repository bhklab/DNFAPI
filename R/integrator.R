## Description

#Our standard network contains 365 drugs in it, with 8 of those drugs having to be
#imputed in the perturbation layer since we're using the new perturbation signatures.
#The aim of this notebook is to create the largest network possible by using all information
#from the sensitivity layer, and not subsetting to just the drugs that are in common with
#structure and perturbation.

integrator <- function(struct, perturb, sens){
  #Get all necessary imports and functions
  set.seed(9833)

  badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

  #Load LINCS metadata file. Used it contains the smiles string for the drugs in the perturbation layer.
  #These are then used to construct the structure layer.

  lincs.meta <- read.csv("../Data/LINCS.csv", stringsAsFactors = FALSE)
  lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
  lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

  #Load in senstivity, imaging data as well as perturbation signatures

  pert.file.name <- "../Data/pert_features_full.RData"
  sensitivity.file.name <- "../Data/combined_sensitivity/combined_sens_iname_replaced.RData"

  sens.data <- SensitivityDataFlexible(sensitivity.file.name)  ## 645 X 239
  pert.data <- PerturbationDataFlexibleCustomSig(pert.file.name)  ## 978 X 23

  #Subset the lincs metadata file to only the drugs appearing in the sensitivity layer

  sens.names <- rownames(sens.data)

  lincs.meta.subset <- lincs.meta[match(sens.names, lincs.meta$pert_iname),]
  lincs.meta.subset <- lincs.meta.subset[!is.na(lincs.meta.subset$X),]

  #Now we create the structure fingerprints based on the subsetted lincs metadata

  strc.data <- StructureDataFlexible(lincs.meta.subset)  ## a vector  --> 239 elemnts

  ## Find Common Drugs

  #Now to find the drugs that are common between all layers

  #Now we can choose the layers to be used in the integration and see the drugs that are common between them.
  #Leave out pert names from the intersection since we know that the pert layer is missing 8 drugs compared to
  #the strc layer so well use our imputation method to get a full integrated network of 365 drugs instead of 357 drugs.

  layers <- list(sens.names = sort(colnames(sens.data)),
                 strc.names = names(strc.data))
  common.drugs <- Reduce(intersect, Filter(Negate(is.null),layers))
  print(length(common.drugs))

  ## Create Similarities

  #At this point we have all the data necessary to create the various similarity matrices. Keep in mind that in the case of the senstivity layer, the data loaded in is already a similarity matrix, so there is nothing to do there. Note that we are leaving out the perturbation layer since that is the one that we are evaluating the effect of imputation on.

  strc.cor <- fingerprint::fp.sim.matrix(strc.data, method = "tanimoto")
  colnames(strc.cor) <- names(strc.data)
  rownames(strc.cor) <- names(strc.data)

  strc.cor <- strc.cor[common.drugs, common.drugs]

  pert.cor <- cor(pert.data[, intersect(common.drugs, colnames(pert.data))], method="pearson", use="pairwise.complete.obs")

  pert.cor <- medianSimilarity(list(pert.cor))[[1]]

  sens.cor <- sens.data

  ## Integrate Network

  all.drugs <- sort(rownames(sens.data))

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
