#' Database searching by MSGF+
#'
#' This part contains MSGF+ database searching modules integrated in UHR-IonStar application.
#'
#' @param path the location containing .raw files and .fasta database
#' @param msgfpar specified parameters for MSGF+
#' @param msgfpath the location containing MSGF+ java file
#' @return identification information .mzid file
#' @export


DBsearching <- function (path, msgfpar, msgfpath){
  withProgress(message = 'Protein Identification', value = 0, {
    incProgress(0.2, detail = '.raw to .mzxml')
    convertMSFiles(
      files = path,
      outPath = NULL,
      dirs = TRUE,
      anaInfo = NULL,
      from = "thermo",
      to = "mzXML",
      overWrite = FALSE,
      algorithm = "pwiz",
      centroid = 'vendor',
      filters = NULL,
      extraOpts = NULL,
      PWizBatchSize = 1
    )
    incProgress(0.4, detail = 'MSGF+ database searching')
    setwd(path)
    files <- list.files(path = path)
    runMSGF(msgfpar, files[grep('.mzXML',files)], memory = 10000, msgfPath = msgfpath)
  })
}


modi_names <- c('Carbamidomethyl','Oxidation','Deamidated','Methyl','Acetyl','Phospho',
                'Glu->pryo-Glu','Gln->pryo-Glu')
modification_list <- list(
  list(name='Carbamidomethyl',
       composition='C2H3N1O1',
       residues='C',
       type='fix',
       position='any'),
  list(name='Oxidation',
       composition='O1',
       residues='M',
       type='opt',
       position='any'),
  list(name='Deamidated',
       composition='H-1N-1O1',
       residues='NQ',
       type='opt',
       position='any'),
  list(name='Methyl',
       composition='CH2',
       residues='K',
       type='opt',
       position='any'),
  list(name='Acetyl',
       composition='C2H2O',
       residues='*',
       type='opt',
       position='Prot-N-term'),
  list(name='Phospho',
       composition='HO3P',
       residues='STY',
       type='opt',
       position='any'),
  list(name='Glu->pryo-Glu',
       composition='H-2O-1',
       residues='E',
       type='opt',
       position='N-term'),
  list(name='Gln->pryo-Glu',
       composition='H-3N-1',
       residues='Q',
       type='opt',
       position='N-term')
)


