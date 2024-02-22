# devtools::install_github("HCBravoLab/metagenomeFeatures")
# install.packages('RSQLite')
# install.packages('DBI')
# install.packages('dplyr')

library(tidytree)
library(ape)
library(metagenomeFeatures)
library(RSQLite)
library(DBI)
library(dplyr)
library(gtools)

set.seed(59049)

### Set up GreenGenes DB ###
gg85 <- metagenomeFeatures::get_gg13.8_85MgDb()
tree <- gg85@tree
metadata <- list(ACCESSION_DATE = date(),
                 URL = "ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus",
                 DB_TYPE_NAME = "GreenGenes",
                 DB_VERSION = "13.8 85% OTUS",
                 DB_TYPE_VALUE = "MgDb",
                 DB_SCHEMA_VERSION = "2.0")

gg_db_file <- system.file("extdata", 'gg13.8_85.sqlite',
                          package = "metagenomeFeatures")

db_conn <- DBI::dbConnect(RSQLite::SQLite(), gg_db_file)

taxa_dbi <- dplyr::tbl(src = db_conn, from = "Seqs")

### Collect subtree entries from DB ###

keys <- "Gammaproteobacteria"
keytype <- "Class"
columns <- "all"

if (!is.null(keys)) {
  if (keytype !=  "Keys") {
    level_id <- stringr::str_sub(string = keytype,
                                 start = 1,
                                 end = 1) %>%
      tolower() %>% rep(length(keys))
    if (metadata$DB_TYPE_NAME == "GreenGenes") {
      keys <- stringr::str_c(level_id,keys,sep = "__")
    }
  }
  select_tbl <- dplyr::filter(taxa_dbi, Class ==  keys)
}else{
  select_tbl <- taxa_dbi
}

if (columns[1] != "all") {
  select_tbl <- dplyr::select(select_tbl, .dots = columns)
}

taxa_df <- dplyr::collect(select_tbl)

### Create subtree ###

ids=taxa_df$Keys
drop_tips <- tree$tip.label[!(tree$tip.label %in% ids)]
tree<-ape::drop.tip(tree,drop_tips) %>% ape::as.phylo()

### Model effect ###

Clade1=tidytree::offspring(tree, 478)
Clade1=Clade1[Clade1<=247]
Clade2=tidytree::offspring(tree, 388)
Clade2=Clade2[Clade2<=247]

f_s <- function(x){sum(sin(2*pi*x))}
f_m <- function(x){mean(cos(2*pi*x))}
f_y <- function(y){cos(pi*y/3)}

condgen <- function(X) {
  f_s <- function(x){sum(sin(2*pi*x))}
  
  Y<-rnorm(1, mean = f_s(X))
  
  Y
}

n <- 100
m <- 5000
beta <- 100
model <- 1

X=matrix(runif(n*5), nrow = n)
means=apply(X,1,f_s)
Y=rnorm(n,mean = means)
P=matrix(0,nrow = length(tree$tip.label)+tree$Nnode,ncol = n)
if (model==1) {
  means=apply(X,1,f_m)*0+f_y(Y)*1
} else if (model==2) {
  means=apply(X,1,f_m)*0.2+f_y(Y)*0.8
} else if (model==3) {
  means=apply(X,1,f_m)*0.4+f_y(Y)*0.6
} else if (model==4) {
  means=apply(X,1,f_m)*0.6+f_y(Y)*0.4
} else if (model==5) {
  means=apply(X,1,f_m)*0.8+f_y(Y)*0.2
}
for (i in 1:n) {
  alphaQ=c(rep(1,length(tree$tip.label)),rep(0,tree$Nnode))
  alphaQ[Clade1]=alphaQ[Clade1]+means[i]
  alphaQ[Clade2]=alphaQ[Clade2]-means[i]
  P[,i] = gtools::rdirichlet(1, alphaQ)
  P[,i] = rmultinom(1,m,P[,i])
  P[,i] = P[,i]/sum(P[,i])
}

ExX=matrix(runif(n*5*beta), nrow = n*beta)
means=apply(ExX,1,f_s)
ExY=rnorm(n*beta,mean = means)
Exk <- 30

### Save data ###

test_ConIndDAFOT <- list(P = P, Y = Y, X = X, tree = tree, condgen = condgen, ExY = ExY, ExX = ExX, Exk = Exk)
usethis::use_data(test_ConIndDAFOT, overwrite = TRUE)
