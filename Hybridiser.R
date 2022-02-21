# Hybridiser v0.1 
# Elisa Somenzi & Mario Barbato 2020
# If you  use this script please cite:
# Identification of Ancestry Informative Marker panels to assess hybridisation between feral and domestic sheep 
# (Somenzi E, Ajmone Marsan P, Barbato M, 2020) DOI: 10.3390/ani10040582
# 
# Parameters:
# parent_1.ped  : ped file for parent population 1
# parent_1.map  : ped file for parent population 1
# parent_2.ped  : ped file for parent population 1
# parent_2.map  : ped file for parent population 1
# number_HYB    : number of hybrid individuals to simulate 
# out_FID       : Family ID to assign in the hybrid .ped file 
# outfile       : file name of the otput files
# miss          : character used for missing data
# 
# Usage       : Hybridiser(pop_1.ped, pop_1.map, pop_2.ped, pop_2.map) 

Hybridiser <- function(
  parent_1.ped = "",
  parent_1.map = "",
  parent_2.ped = "",
  parent_2.map = "",
  number_HYB = 10,
  out_FID = "HYB",
  outfile = "HYB",
  miss = '0'
){
  cat("Hybridiser v0.1\n\n")
  
  require(data.table)
  # check input parameters ----
  # files exist
  if(!file.exists(parent_1.ped)) stop("Could not open parent_1.ped", call. = F)
  if(!file.exists(parent_1.map)) stop("Could not open parent_1.map", call. = F)
  if(!file.exists(parent_2.ped)) stop("Could not open parent_2.ped", call. = F)
  if(!file.exists(parent_2.map)) stop("Could not open parent_2.map", call. = F)
  
  # number hyb > 0
  if(number_HYB < 1) stop("Could not generate less than one hybrid", call. = F)
  
  # out_FAM len > 0
  if(nchar(out_FAM) < 1) stop("Need a non-zero length Family ID", call. = F)
  
  # outfile len > 0
  if(nchar(outfile) < 1) stop("Need a non-zero length output file name", call. = F)
  
  # load files ----
  p1.p <- fread(parent_1.ped, drop = 1:6)
  p1.m <- fread(parent_1.map)
  p2.p <- fread(parent_2.ped, drop = 1:6)
  p2.m <- fread(parent_2.map)
  
  num_loci = NROW(p1.m)
  num_ID1 = NROW(p1.p)
  num_ID2 = NROW(p2.p)
  
  if(num_loci != NROW(p2.m)) stop('the two files need the same number of loci', call. = F)
  
  refAllp1 = character(num_loci)
  refAllp2 = character(num_loci)
  
  altAll1 = character(num_loci)
  altAll2 = character(num_loci)
  
  frqAllp1 = numeric(num_loci)
  frqAllp2 = numeric(num_loci)
  
  c = 0 
  
  for(x in seq(1, NCOL(p1.p) - 1, by = 2)){
    c = c + 1
    
    # ped 1
    listAll <- c(p1.p[[x]], p1.p[[x + 1]])
    t1 = table(listAll[!listAll %in% miss])
    refAllp1[c] = names(t1[1])
    frqAllp1[c] = t1[[1]]/sum(t1)
    altAll1[c] = ifelse(NROW(t1) == 2, names(t1[2]), refAllp1[c])
    
    # ped 2
    listAll <- c(p2.p[[x]], p2.p[[x + 1]])
    t2 = table(listAll[!listAll %in% miss])
    refAllp2[c] = names(t2[1])
    frqAllp2[c] = t2[[1]]/sum(t2)
    altAll2[c] = ifelse(NROW(t2) == 2, names(t2[2]), refAllp2[c])
  }
  
  # Hybridise ----
  df_list <- list()
  
  for(x in 1:number_HYB){
    
    hyb_snp = character(num_loci)
    
    for(snp in 1:num_loci){
      
      a1 = ifelse(runif(1) <= frqAllp1[snp], refAllp1[snp], altAll1[snp])
      a2 = ifelse(runif(1) <= frqAllp2[snp], refAllp2[snp], altAll2[snp])
      
      hyb_snp[snp] = paste(a1, a2)
      
      df_list[[x]] <- paste0(out_FID, ' ', out_FID, '_', x, ' 0 0 0 -9 ',  paste0( hyb_snp, collapse = ' '))
    }
  }
  
  # print ped
  write(unlist(df_list), paste0(outfile, '.ped'))
  # print map
  fwrite(p1.m, paste0(outfile, '.map'), col.names = F, sep = '\t')
  
  cat("Done -", number_HYB , "hybrid IDs generated\n")
  cat("Results copied as ", outfile,".{ped,map}\n", sep = "")
}

