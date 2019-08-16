#' Returns the default integrated matrix as well as impacts of the original three modalities on the integrated matrix
#' @importFrom apcluster apcluster
#'
#' @return The default integrated matrix
#'
#' @export
getPremadeNetwork <- function() {
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


  # Generate the profile for all drugs
  drugNames <- rownames(result)
  profiles <- lapply(drugNames, getDrugInfo)
  names(profiles) <- drugNames

  return(list(impact=impact, result=result, profiles=profiles))
}

#' Returns a network that integrates all data modalities
#' @importFrom stringr str_remove_all
#' @importFrom apcluster apcluster
#' @import SNFtool
#'
#' @return The new network
#'
#' @export
getNewNetwork <- function(...) {
  impact <- list(item=0) # A list containing the influence of each data type on each connection

  # Integrate the data
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

  # Generate the profile for all drugs
  drugNames <- rownames(integrated)
  profiles <- lapply(drugNames, getDrugInfo)

  return(list(impact=impact, result=integrated, profiles=profiles))
}

convertToURL <- function(id, drugName, database){
  if (database == 'CLUE.IO'){
    return(paste0("https://clue.io/command?q=", drugName))
  }
  if (database == 'ChEMBL'){
    return(paste0("https://www.ebi.ac.uk/chembl/compound_report_card/", id, "/"))
  }
  if (database == 'DrugBank'){
    return(paste0("https://www.drugbank.ca/drugs/", id))
  }
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

  #Get all rows from the drug-info dataframe that match the drug name
  allDrugInfo <- drugTargetInfo[which(drugTargetInfo[,'MOLECULE_NAME'] == drug & drugTargetInfo[,'DATABASE'] != 'CTRPv2'), ]
  if (nrow(allDrugInfo) == 0){
    return(NULL)
  }

  drugTargets <- allDrugInfo[['TARGET_NAME']]

  #Get all the links to drug databases
  ids <- allDrugInfo[['ID']]
  links <- mapply(convertToURL, id = ids, drugName=allDrugInfo[['MOLECULE_NAME']], database=allDrugInfo[['DATABASE']])

  return(list(targets=sort(unique(drugTargets)), links=unique(links)))
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
  uniqueDrugs <- gsub(badchars, "", uniqueDrugs)

  return(toupper(uniqueDrugs))
}

#' Test function
#'
#' @return A list
#'
#' @export
# returnList <- function(){
#   jaccuzi <- "split"
#   return(list(list(zerg='zerg', shit=c(12,3,5)), c(1,2,3,4,5,6), nig="yes", $jaccuzi = 'yo'))
# } #I should use names(lst) <- c("x", "y", "z") to set the names of all drug info objects
