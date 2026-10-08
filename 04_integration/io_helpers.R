###############################################################################
# io_helpers.R                                                       2026-10-08
# metilene_input(dir, id): the metilene output for one comparison, plain or
# xz-compressed, in a form fread() accepts (a path, or a decompression command).
###############################################################################
metilene_input <- function(dir, id) {
  f <- file.path(dir, "metilene", paste0(id, "_metilene_DMRs_147_cpn3.tsv"))
  if (file.exists(f)) return(list(file = f))
  if (file.exists(paste0(f, ".xz"))) return(list(cmd = sprintf("xz -dc '%s.xz'", f)))
  stop("no metilene output for ", id, " in ", dir)
}
fread_metilene <- function(dir, id, ...) do.call(data.table::fread, c(metilene_input(dir, id), list(...)))
