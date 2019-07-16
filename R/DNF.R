#' Returns the default integrated matrix as well as impacts of the original three modalities on the integrated matrix
#' @importFrom apcluster apcluster
#'
#' @return The default integrated matrix
#'
#' @export
getBigNetwork <- function() {
  result <- as.data.frame(communityAugment(integrated80, GMT_TARG))
  impact <- list(item=0) # A list containing the influence of each data type on each connection
  # numDrugs <- ncol(result)
  # impact1 <- matrix(rexp(numDrugs * numDrugs, rate=.1), ncol=numDrugs)
  # impact2 <- matrix(rexp(numDrugs * numDrugs, rate=.1), ncol=numDrugs)
  # impact3 <- matrix(rexp(numDrugs * numDrugs, rate=.1), ncol=numDrugs)
  # sum <- impact1 + impact2 + impact3
  # impact1 <- as.data.frame(impact1 / sum)
  # impact2 <- as.data.frame(impact2 / sum)
  # impact3 <- as.data.frame(impact3 / sum)
  # colnames(impact1) <- colnames(result)
  # rownames(impact1) <- colnames(result)
  # colnames(impact2) <- colnames(result)
  # rownames(impact2) <- colnames(result)
  # colnames(impact3) <- colnames(result)
  # rownames(impact3) <- colnames(result)
  return(list(impact=impact, result=result))
}

#' Returns a network that integrates all data modalities
#' @importFrom stringr str_remove_all
#' @importFrom apcluster apcluster
#'
#' @return The new network
#'
#' @export
getNewNetwork <- function(...) {
  dataNames <- names(list(...))
  sensitivity <- as.data.frame(list(...)[dataNames[2]], stringsAsFactors=FALSE)
  colnames(sensitivity) <- str_remove_all(colnames(sensitivity), paste(dataNames[2], ".", sep=""))
  rownames(sensitivity) <- sensitivity$RowName
  sensitivity$RowName <- NULL
  sensitivity <- as.matrix(as.numeric(as.character(sensitivity)))

  perturbation <- as.data.frame(list(...)[dataNames[1]], stringsAsFactors=FALSE)
  colnames(perturbation) <- str_remove_all(colnames(perturbation), paste(dataNames[1], ".", sep=""))
  rownames(perturbation) <- perturbation$RowName
  perturbation$RowName <- NULL
  perturbation <- as.numeric(as.character(perturbation))

  structure <- as.data.frame(list(...)[dataNames[3]], stringsAsFactors=FALSE)
  colnames(structure) <- str_remove_all(colnames(structure), paste(dataNames[3], ".", sep=""))
  rownames(structure) <- structure$RowName
  structure$RowName <- NULL
  integrated <- integrator(structure, perturbation, sensitivity)
  integrated <- as.data.frame(communityAugment(integrated, GMT_TARG))
  impact <- list(item=0) # A list containing the influence of each data type on each connection
  return(list(impact=impact, result=integrated))
}

#' Returns information related to a drug
#'
#' @param drug The name of the drug
#' @return The CHEMBL ID of the drug as well as a list of drug targets
#'
#' @export
getDrugInfo <- function(drug){
  badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  drugTargetInfo$MOLECULE_NAME <- gsub(badchars, "", drugTargetInfo$MOLECULE_NAME)

  drug <- tolower(drug) #turn the drug name to lowercase for easy matching

  #Get all the targets from the dataframe
  allDrugInfo <- drugTargetInfo[which(drugTargetInfo[,'MOLECULE_NAME'] == drug), ]
  if (nrow(allDrugInfo) == 0){
    return(c())
  }

  drugTargets <- allDrugInfo[['TARGET_NAME']]

  #Get the CHEMBL ID
  id <- allDrugInfo[['ID']][1]

  return(append(id, drugTargets))
}

#' Returns the drugs with input targets
#'
#' @param targets The names of the drug targets
#' @return A list of drugs
#'
#' @export
getDrugsOfTargets <- function(targets){

  #drugTargetInfo$MOLECULE_NAME <- gsub(badchars, "", drugTargetInfo$MOLECULE_NAME)

  targets <- toupper(targets) #turn the target name to uppercase for easy matching

  #Get all the targets from the dataframe
  allDrugInfo <- drugTargetInfo[which(drugTargetInfo[,'TARGET_NAME'] %in% targets), ]
  if (nrow(allDrugInfo) == 0){
    return(c())
  }

  uniqueDrugs <- unique(allDrugInfo[['MOLECULE_NAME']])

  badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  uniqueDrugsE <- gsub(badchars, "", uniqueDrugs)

  return(toupper(uniqueDrugs))
}
