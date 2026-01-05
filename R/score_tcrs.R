#' Score TCR Sequences based on Pre-trained Models
#'
#' This function takes a data frame of TCR sequences and calculates scores
#' associated with different T cell phenotypes (e.g., Innate-like, CD8, 
#' Regulatory, Memory). It first featurizes the TCRs and then applies a 
#' weighted scoring model.
#'
#' @param data A data frame containing TCR sequence information. It must contain
#' columns corresponding to the chain(s) being scored. For alpha-beta pairs 
#' (`chain = "ab"`), it requires: `TCRA_vgene`, `TCRA_jgene`, `TCRA_cdr3aa`, 
#' `TCRB_vgene`, `TCRB_jgene`, and `TCRB_cdr3aa`. For single chains, it 
#' requires the respective `TCRA` or `TCRB` columns.
#' @param chain A character string specifying the TCR chain type.
#'   Accepted values are:
#'   \itemize{
#'     \item `"ab"` for paired alpha-beta chains.
#'     \item `"a"` for alpha chain only.
#'     \item `"b"` for beta chain only.
#'   }
#' @param MAIT_NKT A logical value (`TRUE` or `FALSE`). If `TRUE` and `chain = "b"`,
#' the function will apply models to score for MAIT-like and NKT-like phenotypes.
#' Defaults to `FALSE`.
#' 
#' @examples
#' # Score paired alpha-beta TCRs
#' ab_scores <- score_tcrs(tcrpheno_clones[1:200,], 
#'                         chain = "ab")
#' 
#' # Score beta-chain only for MAIT/NKT phenotypes
#' nkt_scores <- score_tcrs(tcrpheno_clones[1:200,], 
#'                          chain = "b", 
#'                          MAIT_NKT = TRUE)
#'
#' @return A data frame containing the calculated scores for each TCR. The row names
#' will correspond to the input TCR identifiers, and column names will indicate the
#' score type (e.g., "TCR-innate", "TCR-CD8").
#' @export
score_tcrs <- function(data, chain, MAIT_NKT = FALSE){
  # Define score names based on chain type
  score_names = c("TCR-innate", "TCR-CD8", "TCR-reg", "TCR-mem")
  ftz = featurize_tcrs(data, chain)
  print("TCRs featurized!")
  # Select appropriate weights and scaling parameters based on the chain type
  if (chain=="ab"){
    weights = weightsAB
    score_means = ABscore_mns
    score_sds = ABscore_sds
  } else if (chain=="a"){
    weights = weightsA
    score_names = gsub("TCR", "TCRalpha", score_names)
    score_means = Ascore_mns
    score_sds = Ascore_sds
  } else if (chain=="b"){
    weights = weightsB
    score_names = gsub("TCR", "TCRbeta", score_names)
    if (MAIT_NKT == TRUE){
      weights = weightsMAITNKT
      score_names = c("TCRbeta-MAIT", "TCRbeta-NKT")
      score_means = MNKTscore_mns
      score_sds = MNKTscore_sds
    } else {
      score_means = Bscore_mns
      score_sds = Bscore_sds
    }
  } else {
    print("please specify the 'chain' argument (a, b, or ab)")
  }
  # Get mean and standard deviation for feature scaling
  m = mns[as.character(rownames(weights)),]
  s = sds[as.character(rownames(weights)),]
  print("scoring TCRs...")
  rownames(ftz) = as.character(ftz$id)
  ftz = scale_variables(ftz, m, s)
  # Calculate scores via matrix multiplication
  scores = as.matrix(ftz) %*% as.matrix(weights)
  rownames(scores) = rownames(ftz)
  colnames(scores) = score_names
  # Set up for final score scaling
  names(score_means) = score_names
  names(score_sds) = score_names
  # Scale the final scores
  scores = scale_variables(scores, score_means, score_sds)
  print("all done!")
  return(data.frame(scores))
}
