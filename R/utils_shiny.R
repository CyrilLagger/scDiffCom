reduce_go_terms <- function(
  object,
  method,
  threshold
) {
  ASPECT <- OR_UP <- BH_P_VALUE_UP <- VALUE_BIS <-
    OR_DOWN <- BH_P_VALUE_DOWN <- OR_FLAT <- BH_P_VALUE_FLAT <- NULL
  message(
    "Reducing GO terms based on semantic similarity"
  )
  if (object@parameters$LRI_species == "human") {
    orgdb <- "org.Hs.eg.db"
  } else  if (object@parameters$LRI_species == "mouse") {
    orgdb <- "org.Mm.eg.db"
  } else  if (object@parameters$LRI_species == "rat") {
    orgdb <- "org.Rn.eg.db"
  } else {
    stop(
      "Species must be 'human', 'mouse' or 'rat'"
    )
  }
  if (!("GO_TERMS" %in% names(object@ora_table))) {
    stop(
      "Table 'GO_TERMS' must be present in slot 'ora_table'"
    )
  }
  semd_BP <- GOSemSim::godata(
    OrgDb = orgdb,
    ont = "BP"
  )
  semd_MF <- GOSemSim::godata(
    OrgDb = orgdb,
    ont = "MF"
  )
  semd_CC <- GOSemSim::godata(
    OrgDb = orgdb,
    ont = "CC"
  )
  all_org_go <- unique(
    c(
      names(semd_BP@IC),
      names(semd_CC@IC),
      names(semd_MF@IC))
  )
  dt <- copy(object@ora_table$GO_TERMS)
  dt <- dt[VALUE_BIS %in% all_org_go]
  rbindlist(
    l = lapply(
      stats::setNames(
        sort(unique(dt$ASPECT)),
        sort(unique(dt$ASPECT))
      ),
      function(aspect) {
        rbindlist(
          l = lapply(
            list(UP = "UP", DOWN = "DOWN", FLAT = "FLAT"),
            function(regulation) {
              dt_intern <- dt[ASPECT == aspect]
              if (regulation == "UP") {
                dt_intern <- dt_intern[OR_UP >= 1 & BH_P_VALUE_UP <= 0.05]
                OR_intern <- dt_intern$OR_UP
                BH_intern <- dt_intern$BH_P_VALUE_UP
              } else if (regulation == "DOWN") {
                dt_intern <- dt_intern[OR_DOWN >= 1 & BH_P_VALUE_DOWN <= 0.05]
                OR_intern <- dt_intern$OR_DOWN
                BH_intern <- dt_intern$BH_P_VALUE_DOWN
              } else if (regulation == "FLAT") {
                dt_intern <- dt_intern[OR_FLAT >= 1 & BH_P_VALUE_FLAT <= 0.05]
                OR_intern <- dt_intern$OR_FLAT
                BH_intern <- dt_intern$BH_P_VALUE_FLAT
              }
              if (nrow(dt_intern) == 0) return(NULL)
              if (any(is.infinite(OR_intern))) {
                if (all(is.infinite(OR_intern))) {
                  OR_intern <- rep(2, length(OR_intern))

                } else {
                  max_finite <- max(OR_intern[is.finite(OR_intern)])
                  OR_intern[is.infinite(OR_intern)] <- max_finite
                }
              }
              scores_intern <- -log10(BH_intern) * log2(OR_intern)
              scores_intern <- stats::setNames(
                scores_intern,
                dt_intern$VALUE_BIS
              )
              if (aspect == "biological_process") {
                ont_intern <- "BP"
                semd <- semd_BP
              } else if (aspect == "molecular_function") {
                ont_intern <- "MF"
                semd <- semd_MF
              } else if (aspect == "cellular_component") {
                ont_intern <- "CC"
                semd <- semd_CC
              }
              simMatrix <- rrvgo::calculateSimMatrix(
                x = dt_intern$VALUE_BIS,
                orgdb = orgdb,
                semdata = semd,
                ont = ont_intern,
                method = method
              )
              if (is.null(dim(simMatrix))) return(NULL)
              reducedTerms <- rrvgo::reduceSimMatrix(
                simMatrix = simMatrix,
                scores = scores_intern,
                threshold = threshold,
                orgdb = orgdb
              )
              if (nrow(reducedTerms) == 0) return(NULL)
              setDT(reducedTerms)
              reducedTerms
            }
          ),
          idcol = "REGULATION"
        )
      }
    ),
    idcol = "ASPECT"
  )
}

validate_reduced_go_table <- function(
  object,
  reduced_go_table
) {
  REGULATION <- VALUE_BIS <- OR_UP <- BH_P_VALUE_UP <-
    OR_DOWN <- BH_P_VALUE_DOWN <- OR_FLAT <- BH_P_VALUE_FLAT <- NULL
  orig_go_dt <- copy(object@ora_table$GO_TERMS)
  check_1 <- all(reduced_go_table[REGULATION == "UP"]$go %in%
                   orig_go_dt[OR_UP >=1 & BH_P_VALUE_UP <= 0.05]$VALUE_BIS) &
    all(reduced_go_table[REGULATION == "DOWN"]$go %in%
          orig_go_dt[OR_DOWN >=1 & BH_P_VALUE_DOWN <= 0.05]$VALUE_BIS) &
    all(reduced_go_table[REGULATION == "FLAT"]$go %in%
          orig_go_dt[OR_FLAT >=1 & BH_P_VALUE_FLAT <= 0.05]$VALUE_BIS)
  if(!check_1) {
    return(FALSE)
  }
  dt_UP <- orig_go_dt[
    OR_UP >= 1 &
      BH_P_VALUE_UP <= 0.05 &
      VALUE_BIS %in% reduced_go_table[REGULATION == "UP"]$go
  ]
  OR_UP_intern <- dt_UP$OR_UP
  dt_DOWN <- orig_go_dt[
    OR_DOWN >=1 &
      BH_P_VALUE_DOWN <= 0.05 &
      VALUE_BIS %in% reduced_go_table[REGULATION == "DOWN"]$go
  ]
  OR_DOWN_intern <- dt_DOWN$OR_DOWN
  dt_FLAT <- orig_go_dt[
    OR_FLAT >= 1 &
      BH_P_VALUE_FLAT <= 0.05 &
      VALUE_BIS %in% reduced_go_table[REGULATION == "FLAT"]$go
  ]
  OR_FLAT_intern <- dt_FLAT$OR_FLAT
  if (any(is.infinite(OR_UP_intern))) {
    if (all(is.infinite(OR_UP_intern))) {
      OR_UP_intern <- rep(2, length(OR_UP_intern))
    } else {
      max_finite <- max(OR_UP_intern[is.finite(OR_UP_intern)])
      OR_UP_intern[is.infinite(OR_UP_intern)] <- max_finite
    }
  }
  if (any(is.infinite(OR_DOWN_intern))) {
    if (all(is.infinite(OR_DOWN_intern))) {
      OR_DOWN_intern <- rep(2, length(OR_DOWN_intern))
    } else {
      max_finite <- max(OR_DOWN_intern[is.finite(OR_DOWN_intern)])
      OR_DOWN_intern[is.infinite(OR_DOWN_intern)] <- max_finite
    }
  }
  if (any(is.infinite(OR_FLAT_intern))) {
    if (all(is.infinite(OR_FLAT_intern))) {
      OR_FLAT_intern <- rep(2, length(OR_FLAT_intern))
    } else {
      max_finite <- max(OR_FLAT_intern[is.finite(OR_FLAT_intern)])
      OR_FLAT_intern[is.infinite(OR_FLAT_intern)] <- max_finite
    }
  }
  check_2_UP <- identical(
    sort(-log2(OR_UP_intern)*log10(dt_UP$BH_P_VALUE_UP)),
    sort(reduced_go_table[REGULATION == "UP"]$score)
  )
  check_2_DOWN <- identical(
    sort(-log2(OR_DOWN_intern)*log10(dt_DOWN$BH_P_VALUE_DOWN)),
    sort(reduced_go_table[REGULATION == "DOWN"]$score)
  )
  check_2_FLAT <- identical(
    sort(-log2(OR_FLAT_intern)*log10(dt_FLAT$BH_P_VALUE_FLAT)),
    sort(reduced_go_table[REGULATION == "FLAT"]$score)
  )
  if (!all(check_2_UP, check_2_DOWN, check_2_FLAT)) {
    return(FALSE)
  }
  TRUE
}

