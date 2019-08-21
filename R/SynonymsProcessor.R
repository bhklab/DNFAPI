getUnique <- function(x){
  aliases <- as.character(drugs_with_ids[x, ])
  aliases <- aliases[!is.na(aliases)]
  aliases <- aliases[!grepl('/', aliases)]
  aliases <- toupper(aliases)
  return(unique(aliases))
}
