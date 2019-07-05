getSmallNetwork <- function() {
  result <- applyCommunity(default_integrated_matrix)[0:100, 0:100]
  impact <- list(item=0) # A list containing the influence of each data type on each connection
  return(list(impact=impact, result=result))
}

#' Returns the default integrated matrix as well as impacts of the original three modalities on the integrated matrix
#'
#' @importFrom grDevices rgb2hsv
#'
#' @return The default integrated matrix
#'
#' @export
getBigNetwork <- function() {
  result <- applyCommunity(default_integrated_matrix)
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

#' Returns the a small network
#'
#' @param newData User's own data in JSON format (converted to an R data type)
#'
#' @return The new network
#'
#' @export
getNewNetwork <- function(newData) {
  return(getSmallNetwork())
}

jsonInput <- function(input){
  return(input[1,1])
}

#' Test
#' @import stringr
#'
#' @param ... A list of JSON-formatted objects.
#'
#' @return Test
#'
#' @export
multiJsonInput <- function(...){
  dataNames <- names(list(...))
  structure <- as.data.frame(list(...)[dataNames[1]])
  perturbation <- as.data.frame(list(...)[dataNames[2]])
  colnames(perturbation) <- str_remove_all(colnames(perturbation), paste(dataNames[2], ".", sep=""))
  sensitivity <- as.data.frame(list(...)[dataNames[3]])
  return(list(structure, perturbation, sensitivity))
}

readFile <- function(name){
  return(read.csv(name))
}

#Add 1 to all drugs that interact with a representative community drug
applyCommunity <- function(data){
  finalDrugs <- c()
  for (row in 1:nrow(data)) {
    community <- as.character(communities[row, 2:ncol(communities)])
    representative <- communities[row,2]
    drugPositions <- which(colnames(default_integrated_matrix) %in% community)
    data[representative, drugPositions] <- data[representative, drugPositions] + 1
    finalDrugs <- append(finalDrugs, drugPositions)
  }
  finalDrugs <- sort(unique(finalDrugs))
  return(data[finalDrugs, finalDrugs])
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
