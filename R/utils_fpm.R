# analyze_FreqItemSets <- function(
#   scdiffcom_dt,
#   target = "closed frequent itemsets",
#   support = 0.05,
#   confidence = 0.1
# ) {
#   transactions <- convert_data_to_transactions(
#     scdiffcom_dt = scdiffcom_dt
#   )
#   res <- compute_freqitemsets_and_rules(
#     transactions = transactions,
#     target = target,
#     support = support,
#     confidence = confidence
#   )
#   freqsets <- res[["freqsets"]]
#   rules <- res[["rules"]]
#   return(list(rules = rules, transactions = transactions))
#   interesting_rules <- get_interesting_rules(
#     rules = rules,
#     transactions = transactions
#   )
#   return(interesting_rules)
# }
#
# convert_data_to_transactions <- function(
#   scdiffcom_dt
# ) {
#   #need to be changed for the correct number of ligand-receptor
#   cols_for_items = c(
#     "L_CELLTYPE", "R_CELLTYPE",
#     "LIGAND_1", "LIGAND_2",
#     "RECEPTOR_1", "RECEPTOR_2", "RECEPTOR_3",
#     "DIFFERENTIALLY_EXPRESSED",
#     "DIFFERENTIAL_DIRECTION"
#   )
#   dt_p <- scdiffcom_dt[, cols_for_items, with = FALSE]
#   dt_p[, DIFFERENTIALLY_EXPRESSED :=
#          fifelse(DIFFERENTIALLY_EXPRESSED, "Sign", "Notsign")]
#   transactions = as(dt_p, "transactions")
#   return(transactions)
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
#       minlen = 1,
#       maxlen = 5,
#       target = target,
#       ext = TRUE,
#       smax = 1,
#       minval = 0,
#       maxtime = 0
#     ),
#     control = list(verbose=FALSE)
#   )
#   rules = arules::ruleInduction(
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
