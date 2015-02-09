### R PACKAGE RCurl REQUIRED ###
rm(list = ls())
library("RCurl")
sys_tst <- .Platform["OS.type"]

### user selects their otu table ###
### THIS MUST BE FORMATED AS TAB DELIMITED TEXT ###
### THE FIRST COLUM SHOULD CONTAIN THE OTU IDS ###
### THE FIRST ROW MUST CONTAIN THE SAMPLE IDS ###
### THUS REMOVE ROWS WITH NOTATIONS GENERATED IN PROGRAMS SUCH AS QIIME E.G., '# Constructed from biom file' ###
### IN THE FIRST ROW, THE COLUMN SHOWING THE TAXONOMIC STRINGS MUST BE LABELED "taxonomy" - ALL LOWERCASE ###
### REMOVE ALL '#' SYMBOLS FROM LABELS - E.G., '#OTU ID' SHOULD BE CHANGED TO 'OTU ID' ###
pth_otu <- file.choose()
otu_table <- read.table(pth_otu, header=T, sep="\t", na.strings="NA")
heads <- names(otu_table)
heads_len <- length(heads)
heads[heads_len+1] <- "function"
heads[heads_len+2] <- "type"
heads[heads_len+3] <- "likelihood"
heads[heads_len+4] <- "notes"
heads[heads_len+5] <- "level"
heads_tlen <- heads_len+5
t_spot <- grep('taxonomy', heads)
tax_strings <- as.vector(otu_table[,t_spot])
Y <- length(tax_strings)
tax_strings <- gsub("_", "@", tax_strings)
tax_strings <- gsub(" ", "@", tax_strings)
tax_strings <- gsub(";", "@", tax_strings)
tax_strings <- gsub(",", "@", tax_strings)
for(i in 1:Y){
  tax_strings[i] <- paste(tax_strings[i], "@", sep = "")
}

if(sys_tst == "windows"){
  ### create the output file path windows ###
  pth_out <- pth_otu
  pth_out <- strsplit(pth_out, split="\\\\")
  pth_out <- pth_out[[1]]
  pth_out_len <- length(pth_out)
  pth_out[pth_out_len] <- "fungal_guild_otu_table.csv"
  pth_out <- paste(pth_out, collapse="\\")
}else{
  ### create the output file path mac ###
  pth_out <- pth_otu
  pth_out <- strsplit(pth_out, split="/")
  pth_out <- pth_out[[1]]
  pth_out_len <- length(pth_out)
  pth_out[pth_out_len] <- "fungal_guild_otu_table.csv"
  pth_out <- paste(pth_out, collapse="/")
}

### read fungal_guild_table.txt from github and parse into a matrix ###
wdata <- getURL("https://raw.githubusercontent.com/UMNFuN/fungal_function/master/fungal_guild_table.txt", .encoding = "UTF-8", ssl.verifypeer = F)
wdata <- strsplit(wdata,"\\n")
wdata <- wdata[[1]]
wdata_len <- length(wdata)
wheads <- c("Taxon", "Taxon Level", "Function", "Type", "Likelihood", "Notes", "Citation/Source")
fdata <- t(matrix(wheads))
for (m in 2:wdata_len){
  m_lin <- wdata[m]
  m_lin <- strsplit(m_lin, "\\t")
  m_lin <- m_lin[[1]]
  fdata <- rbind(fdata, m_lin, deparse.level = 0)
}
ftaxa <- as.vector(fdata[,1])
flevl <- as.vector(fdata[,2])
funct <- as.vector(fdata[,3])
ftype <- as.vector(fdata[,4])
fprob <- as.vector(fdata[,5])
fnote <- as.vector(fdata[,6])
X <- length(ftaxa)

### create output otu table ###
ftab <- matrix(ncol = heads_tlen)

### search for fungal function taxa and build table with discovered taxa ###
for (j in 1:X){
  s_tax <- ftaxa[j]
  s_tax <- gsub(" ", "@", s_tax)
  s_tax <- paste("@", s_tax, "@", sep = "")
  s_spot <- grep(s_tax, tax_strings)
  s_spot_len <- length(s_spot)
  s_fun <- funct[j]
  s_typ <- ftype[j]
  s_pro <- fprob[j]
  s_not <- fnote[j]
  s_lev <- flevl[j]
  inrow <- c()
  if (s_spot_len > 0){
    for (k in 1:s_spot_len){
      inrow <- c()
      r_spot <- s_spot[k]
      for (l in 1:heads_len){
        inrow[l] <- as.character(otu_table[r_spot,l])
      }
      inrow[heads_len+1] <- s_fun
      inrow[heads_len+2] <- s_typ
      inrow[heads_len+3] <- s_pro
      inrow[heads_len+4] <- s_not
      inrow[heads_len+5] <- s_lev
      ftab <- rbind(ftab,inrow)
    }
  }else blnk <- 0
}

### order and dereplicate results table ###
### selecting the most taxonomically resolved duplicates ###
ftab <- ftab[-1,]
ftab2 <- ftab[do.call(order, as.data.frame(ftab)),]
ftab2 <- ftab2[!duplicated(ftab2[,1],fromLast = T),]
ftab2 <- rbind(heads,ftab2)
ftab2 <- ftab2[,(-1*heads_tlen)]
write.table(ftab2, file = pth_out, row.names=F, col.names = F, sep=",")
print("Code Completed")
