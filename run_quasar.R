library(QuASAR)
library(experimentr)
opts <- Experiment.start(list())

ase.dat <- UnionExtractFields(opts$infiles, combine=TRUE)
ase.dat.gt <- PrepForGenotyping(ase.dat = ase.dat, min.coverage=5)

if (length(opts$infiles) == 1) {
  ase.joint <- fitAseNull(ase.dat.gt$ref, ase.dat.gt$alt, log.gmat=log(ase.dat.gt$gmat))

} else {
  ase.joint <- fitAseNullMulti(ase.dat.gt$ref, ase.dat.gt$alt, log.gmat=log(ase.dat.gt$gmat))
}

ourInferenceData <- aseInference(gts=ase.joint$gt, 
                                 eps.vect=ase.joint$eps, 
                                 priors=ase.dat.gt$gmat, 
                                 ref.mat=matrix(ase.dat.gt$ref, ncol = ncol(ase.dat$ref)),
                                 alt.mat=matrix(ase.dat.gt$alt, ncol = ncol(ase.dat$alt)),
                                 min.cov=10, 
                                 sample.names=c(opts$infiles), 
                                 annos=ase.dat.gt$annotations)

for (sample in length(ourInferenceData)) {
  write.table(ourInferenceData[[sample]]$dat,
              file = opts$outfiles[sample],
              quote = FALSE,
              sep="\t",
              row.names=FALSE)
}

Experiment.stop()