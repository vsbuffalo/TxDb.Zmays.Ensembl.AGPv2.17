###
### Load any db objects whenever the package is loaded.
### Copied over from TxDb.Athaliana.BioMart.plantsmart19

.onLoad <- function(libname, pkgname)
{
  ns <- asNamespace(pkgname)
  path <- system.file("extdata", package=pkgname, lib.loc=libname)
  files <- dir(path, pattern="*.sqlite")
  for (i in seq_len(length(files))) {
    db <- loadDb(system.file("extdata", files[[i]], package=pkgname,
                  lib.loc=libname), packageName=pkgname)
    objname <- sub(".sqlite$", "", files[[i]])
    assign(objname, db, envir=ns)
    namespaceExport(ns, objname)
  }
}
