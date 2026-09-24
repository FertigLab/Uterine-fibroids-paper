dom.res <- list.dirs(".")
dom.res <- dom.res[grep("_", dom.res)]

dom.list <- lapply(dom.res, function(res_dir) {
  fname <- list.files(res_dir)[grep("domino.rds", list.files(res_dir))]
  dom <- readRDS(file.path(res_dir, fname))
})

names(dom.list) <- basename(dom.res)

saveRDS(dom.list, paste0("dom_list.cell_type_SnC.rds"))
