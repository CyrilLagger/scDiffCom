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
  ASPECT <- REGULATION <- VALUE_BIS <- OR_UP <- BH_P_VALUE_UP <-
    OR_DOWN <- BH_P_VALUE_DOWN <- OR_FLAT <- BH_P_VALUE_FLAT <- NULL
  orig_go_dt <- copy(object@ora_table$GO_TERMS)
  check_1 <- all(reduced_go_table[REGULATION == "UP"]$go %in%
                   orig_go_dt[OR_UP >= 1 & BH_P_VALUE_UP <= 0.05]$VALUE_BIS) &
    all(reduced_go_table[REGULATION == "DOWN"]$go %in%
          orig_go_dt[OR_DOWN >= 1 & BH_P_VALUE_DOWN <= 0.05]$VALUE_BIS) &
    all(reduced_go_table[REGULATION == "FLAT"]$go %in%
          orig_go_dt[OR_FLAT >= 1 & BH_P_VALUE_FLAT <= 0.05]$VALUE_BIS)
  if(!check_1) {
    return(FALSE)
  }
  orig_go_dt <- orig_go_dt[VALUE_BIS %in% reduced_go_table$go]
  check2 <- lapply(
    sort(unique(orig_go_dt$ASPECT)),
    function(aspect) {
      lapply(
        list(UP = "UP", DOWN = "DOWN", FLAT = "FLAT"),
        function(regulation) {
          dt_intern <- orig_go_dt[ASPECT == aspect]
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
          if (any(is.infinite(OR_intern))) {
            if (all(is.infinite(OR_intern))) {
              OR_intern <- rep(2, length(OR_intern))
            } else {
              max_finite <- max(OR_intern[is.finite(OR_intern)])
              OR_intern[is.infinite(OR_intern)] <- max_finite
            }
          }
          s1 <- sort(
            reduced_go_table[
              REGULATION == regulation & ASPECT == aspect
            ]$score
          )
          if (identical(s1, numeric(0))) return(TRUE)
          identical(
            sort(-log10(BH_intern) * log2(OR_intern)),
            s1
          )
        }
      )
    }
  )
  check2 <- all(unlist(check2))
  if (!check2) return(FALSE)
  TRUE
}

