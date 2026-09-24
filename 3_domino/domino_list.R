dom.res <- list.dirs(".")
dom.res <- dom.res[grep("_", dom.res)]

dom.list <- lapply(dom.res, function(res_dir) {
  fname <- list.files(res_dir)[grep("domino.rds", list.files(res_dir))]
  dom <- readRDS(file.path(res_dir, fname))
})

names(dom.list) <- basename(dom.res)

saveRDS(dom.list, paste0("dom_list.cell_type_SnC.rds"))

domList <- readRDS("dom_list.cell_type_SnC.rds")
sapply(domList, function(dom) {
  table(dom@clusters)
})
meta_df <- data.frame("ID" = c("FP10_Fibroid", "FP10_Myometrium", "FP16_Fibroid", "FP16_Myometrium", "FP6_Fibroid", "FP8_Fibroid", "FP8_Myometrium"), 
                      "group" = c("Fibroid", "Myometrium", "Fibroid", "Myometrium", "Fibroid", "Fibroid", "Myometrium"))

linkage_summary <- summarize_linkages(domList, subject_meta = meta_df)
saveRDS(linkage_summary, "linkage_summary.rds")