# analyze_FreqItemSets <- function(
#   cci_detected,
#   #item_types = c("EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",
#   #               "LIGAND", "RECEPTOR", "REGULATION_SIMPLE",
#   #               "GO_TERMS", "KEGG_PWS"),
#   species,
#   target = "closed frequent itemsets",
#   support = 0.05,
#   confidence = 0.1
# ) {
#   cci_detected <- copy(cci_detected)
#   transactions <- convert_data_to_transactions(
#     cci_detected = cci_detected,
#     species = species
#   )
#   res <- compute_freqitemsets_and_rules(
#     transactions = transactions,
#     target = target,
#     support = support,
#     confidence = confidence
#   )
#   freqsets <- res[["freqsets"]]
#   rules <- res[["rules"]]
#   return(list(rules = rules, freqsets = freqsets, transactions = transactions))
#   interesting_rules <- get_interesting_rules(
#     rules = rules,
#     transactions = transactions
#   )
#   return(interesting_rules)
# }
#
# convert_data_to_transactions <- function(
#   cci_detected,
#   species
# ) {
#   if (species == "mouse") {
#     LRdb_go <- scDiffCom::LRdb_mouse$LRdb_curated_GO
#     LRdb_kegg <- scDiffCom::LRdb_mouse$LRdb_curated_KEGG
#   }
#   if (species == "human") {
#     LRdb_go <- scDiffCom::LRdb_human$LRdb_curated_GO
#     LRdb_kegg <- scDiffCom::LRdb_human$LRdb_curated_KEGG
#   }
#   cci_detected[, LIGAND := paste(LIGAND_1, LIGAND_2, sep =  "_")]
#   cci_detected[, LIGAND := gsub("_NA", "", LIGAND )]
#   cci_detected[, RECEPTOR := paste(RECEPTOR_1, RECEPTOR_2, RECEPTOR_3, sep =  "_")]
#   cci_detected[, RECEPTOR := gsub("_NA", "", RECEPTOR )]
#   tra_list <- lapply(
#     1:nrow(cci_detected),
#     function(i) {
#       temp_go <- LRdb_go[LR_GENES == cci_detected[i]$LR_GENES]$GO_ID
#       temp_kegg <- LRdb_kegg[LR_GENES == cci_detected[i]$LR_GENES]$KEGG_ID
#       res <- c(
#         paste0("EMITTER_CELLTYPE=", cci_detected[i]$EMITTER_CELLTYPE),
#         paste0("RECEIVER_CELLTYPE=", cci_detected[i]$RECEIVER_CELLTYPE),
#         paste0("LIGAND=", cci_detected[i]$LIGAND),
#         paste0("RECEPTOR=", cci_detected[i]$RECEPTOR),
#         paste0("REGULATION_SIMPLE=", cci_detected[i]$REGULATION_SIMPLE)
#       )
#       if (length(temp_go) > 0) {
#         #res <- c(res, paste0("GO=", temp_go))
#       }
#       if (length(temp_kegg) > 0) {
#         #res <- c(res, paste0("KEGG=", temp_kegg))
#       }
#       return(res)
#     }
#   )
#   tra <- as(tra_list, "transactions")
#   temp_iteminfo <- arules::itemInfo(tra)
#   temp_iteminfo$level <- sub("\\=.*", "", temp_iteminfo$labels)
#   temp_iteminfo$value <- sub(".*\\=", "", temp_iteminfo$labels)
#   arules::itemInfo(tra) <- temp_iteminfo
#   #cols_for_items = c(
#   #  "EMITTER_CELLTYPE", "RECEIVER_CELLTYPE",
#   #  "LIGAND", "RECEPTOR",
#   #  "REGULATION_SIMPLE"
#   #)
#   #dt_p <- cci_detected[, cols_for_items, with = FALSE]
#   #transactions = as(dt_p, "transactions")
#   return(tra)
# }
#
# compute_freqitemsets_and_rules <- function(
#   transactions,
#   target,
#   support,
#   confidence
# ) {
#   freqsets <- arules::apriori(
#     data = transactions,
#     parameter = list(
#       support = support,
#       minlen = 2,
#       maxlen = 5,
#       target = target,
#       ext = TRUE,
#       smax = 1,
#       minval = 0,
#       maxtime = 0
#     ),
#     control = list(verbose=FALSE)
#   )
#   rules <- arules::ruleInduction(
#     x = freqsets,
#     transactions = transactions,
#     confidence = confidence,
#     control = list(method = "ptree", verbose = FALSE)
#   )
#   return(list(freqsets = freqsets, rules = rules))
# }
#
# get_interesting_rules <- function(
#   rules,
#   transactions
# ) {
#   sub1 <- get_subset(
#     rules = rules,
#     receptorcell_ligandmolecule=TRUE
#     )
#   sub1_measure = add_interest_measure(sub1, transactions)
#   sub1_dt = convert_to_datatable(sub1, sub1_measure)
#   sub2 = get_subset(rules, ligandcell_receptormolecule=TRUE)
#   sub2_measure = add_interest_measure(sub2, transactions)
#   sub2_dt = convert_to_datatable(sub2, sub2_measure)
#   subsets = list()
#   subsets[["sub1"]] = sub1_dt
#   subsets[["sub2"]] = sub2_dt
#   return(subsets)
# }
#
# get_subset = function(
#   rules,
#   ligandcell_receptormolecule=FALSE,
#   receptorcell_ligandmolecule=FALSE
# ) {
#   if( !xor(ligandcell_receptormolecule, receptorcell_ligandmolecule) ){
#     stop('get_subset: ligandcell_receptormolecule xor receptorcell_ligandmolecule must be true.')
#   }
#   SIGN_DIFF_EXPR_STRING = glue("{s}=Sign", s= "DIFFERENTIALLY_EXPRESSED")
#   # TISSUE_STRING = glue("{s}=", s= "TISSUE")
#   LIGAND1_STRING = glue("{s}=", s= "LIGAND_1")
#   LIGAND2_STRING = glue("{s}=", s= "LIGAND_2")
#   RECEPTOR1_STRING = glue("{s}=", s= "RECEPTOR_1")
#   RECEPTOR2_STRING = glue("{s}=", s= "RECEPTOR_2")
#   RECEPTOR3_STRING = glue("{s}=", s= "RECEPTOR_3")
#   LIGAND_CELLTYPE_STRING = glue("{s}=", s= "L_CELLTYPE")
#   RECEPTOR_CELLTYPE_STRING = glue("{s}=", s= "R_CELLTYPE")
#
#   if (receptorcell_ligandmolecule) {
#     # For ligands - receptor cells
#     sub = subset(
#       rules,
#       subset = rhs %in% SIGN_DIFF_EXPR_STRING
#       # & lhs %pin% TISSUE_STRING
#       & (lhs %pin% LIGAND1_STRING
#          | lhs %pin% LIGAND2_STRING)
#       & lhs %pin% RECEPTOR_CELLTYPE_STRING
#       # & !(lhs %pin% "Direction=")
#     )
#     return(sub)
#   }
#   if (ligandcell_receptormolecule) {
#     # For receptors - ligand cells
#     sub = subset(
#       rules,
#       subset = rhs %in% SIGN_DIFF_EXPR_STRING
#       # & lhs %pin% TISSUE_STRING
#       & (lhs %pin% RECEPTOR1_STRING
#          | lhs %pin% RECEPTOR2_STRING
#          | lhs %pin% RECEPTOR3_STRING)
#       & lhs %pin% LIGAND_CELLTYPE_STRING
#       # & !(lhs %pin% "Direction=")
#     )
#     return(sub)
#   }
# }
#
# add_interest_measure = function(sub, transactions) {
#   measure = interestMeasure(sub,
#                             measure=c("fishersExactTest", "oddsRatio"),
#                             transactions=transactions,
#                             reuse=TRUE)
#   return(measure)
# }
#
# convert_to_datatable = function(sub, sub_measure) {
#   dt = data.table(
#     lhs=labels(lhs(sub)),
#     rhs=labels(rhs(sub)),
#     sub@quality,
#     sub_measure
#   )
#   return(dt)
# }
