## Directory path is indicated with "dirpath" character string
##
## I/O files are in internal R format (RData) using load/save commands
## except "groups.txt", which is used to set width of parallelization.
##
## I/O files are on AFS in
## /u/y/a/yandell/public/html/sysgen/qtlnet/condor
## or equivalently via URL
## http://www.stat.wisc.edu/~yandell/sysgen/qtlnet/condor
##
####################################################################################
parallel.qtlnet <- function(dirpath = "/u/y/a/yandell/public/html/sysgen/qtlnet/condor",
                            phase, ...)
{
  switch(phase,
         qtlnet.phase1(dirpath, ...),
         qtlnet.phase2(dirpath, ...),
         qtlnet.phase3(dirpath, ...),
         qtlnet.phase4(dirpath, ...),
         stop("only 4 phases installed"))
}
####################################################################################
qtlnet.phase1 <- function(dirpath,
                         cross = Pscdbp,
                         pheno.col = 1:13,
                         max.parents = 12,
                         threshold = 3.83,
                         step = 1, ...)
{
  ## PHASE 1: Preparation. Fast. Needed in phases 2 and 3.
  ##
  ## Input files:
  ##       cross.RData
  ## Output files:
  ##       Phase1.RData
  ##       groups.txt
  ##
  ## The "groups.txt" file is used to determine Phase2 width.
  
  ## Calculate genotype probabilities if not already done.
  if (!("prob" %in% names(cross$geno[[1]]))) {
    warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross, step = step)
  }
  
  ## Break up into groups to run on several machines.
  ## 54 groups of ~1000, for a total of 53248 scanone runs.
  parents <- parents.qtlnet(pheno.col, max.parents)
  groups <- group.qtlnet(parents = parents, group.size = 1000)
  
  ## Save all relevant objects for later phases.
  save(cross, pheno.col, max.parents, threshold, parents, groups,
       file = file.path(dirpath, "Phase1.RData"), compress = TRUE)
  
  ## Need to write a file with n.groups lines and group.size columns.
  write.table(groups,
              file = file.path(dirpath, "groups.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  invisible()
}
####################################################################################
qtlnet.phase2 <- function(dirpath, groups.1, groups.2, ...)
{
  ## PHASE 2: Compute BIC scores. Parallelize.
  ##
  ## Input files:
  ##       Phase1.RData
  ##
  ## Output file (one per invocation):
  ##       bicXXX.RData
  ##
  ## The "groups.txt" file (created in Phase1) is used to determine Phase2 width.
  ## That is, groups.txt has 54 lines, hence 54 separate condor runs
  ## to produce bic1.RData through bic54.RData.
  
  ## Assume PERL script gives values groups.1 and groups.2
  ## from a line of the file "groups.txt".

  ## Load Phase 1 computations.
  load(file.path(dirpath, "Phase1.RData"))

  if(missing(groups.1) | missing(groups.2))
    stop("groups.1 and groups.2 must be supplied (from line of file groups.txt")

  ## Pre-compute BIC scores for selected parents.
  bic <- bic.qtlnet(cross, pheno.col, threshold,
                    max.parents = max.parents,
                    parents = parents[seq(groups.1, groups.2)],
                    ...)
  
  save(bic,
       file = paste(tempfile("bic", dirpath), ".RData", sep = ""),
       compress = TRUE)

  invisible()
}
####################################################################################
qtlnet.phase3 <- function(dirpath, ...)
{
  ## PHASE 3: Sample Markov chain (MCMC). Parallelize.
  ##
  ## Input files:
  ##       Phase1.RData
  ##       bicXXX.RData (actual names will use temporary strings in place of XXX)
  ##
  ## Output files (one per invocation):
  ##       mcmcXXX.RData (actual names will use temporary strings for XXX)
  ##
  ## See Phase2 for explanation of "bic*.RData" files.
  ## All "bic.*RData" files in "dirpath" are combined.

  ## Load Phase 1 computations.
  load(file.path(dirpath, "Phase1.RData"))

  ## This could be done once, but it would require splitting this phase in two.
  ## Besides, it is quite fast.
  ## Read in saved BIC scores and combine into one object.
  filenames <- list.files(dirpath, "bic.*RData")
  if(!length(filenames))
    stop(paste("no bic RData files found in", dirpath))
  
  bic.group <- list()
  for(i in seq(length(filenames))) {
    load(file.path(dirpath, filenames[i]))
    bic.group[[i]] <- bic
  }
  saved.scores <- bic.join(cross, pheno.col, bic.group)

  ## This is the phase to parallelize.
  ## Run MCMC with randomized initial network.
  mcmc <- mcmc.qtlnet(cross, pheno.col, threshold = threshold,
                      max.parents = max.parents,
                      saved.scores = saved.scores, init.edges = NULL,
                      ...)
  save(mcmc,
       file = paste(tempfile("mcmc", dirpath), ".RData", sep = ""),
       compress = TRUE)

  invisible()
}
####################################################################################
qtlnet.phase4 <- function(dirpath, ...)
{
  ## PHASE 4: Combine results for post-processing.
  ##
  ## Input files (results of Phase 3):
  ##       mcmcXXX.RData (actual names will use temporary strings in place of XXX
  ##
  ## Output files:
  ##       out.qtlnet.RData
  ##
  ## All "mcmc.*RData" files in "dirpath" are combined.

  ## Combine outputs together.
  filenames <- list.files(dirpath, "mcmc.*RData")
  if(!length(filenames))
    stop(paste("no mcmc RData files found in", dirpath))
  
  outs.qtlnet <- list()
  for(i in seq(length(filenames))) {
    load(file.path(dirpath, filenames[i]))
    outs.qtlnet[[i]] <- mcmc
  }
  
  out.qtlnet <- c.qtlnet(outs.qtlnet)

  save(out.qtlnet,
       file = file.path(dirpath, "out.qtlnet.RData"),
       compress = TRUE)

  out.qtlnet
}
