#' Calculate ROC (Receiver Operating Characteristic) statistics for two numeric vectors.
#'
#' This function calculates ROC statistics, including True Positive Rate (TPR), False Positive Rate (FPR), True Positives (TP), True Negatives (TN), False Positives (FP), False Negatives (FN), and True Skill Score (TSS) for two numeric vectors.
#'
#' @param state A numeric vector representing the state values.
#' @param indicator A numeric vector representing the indicator values.
#' 
#' @return A data frame containing the following columns:
#'   \describe{
#'     \item{state}{The state values.}
#'     \item{label}{A logical vector indicating whether each state is greater than 1 (TRUE) or not (FALSE).}
#'     \item{indicator}{The indicator values.}
#'     \item{TPR}{The True Positive Rate (TPR) calculated as TP / (TP + FN).}
#'     \item{FPR}{The False Positive Rate (FPR) calculated as FP / (FP + TN).}
#'     \item{TP}{The True Positives (TP).}
#'     \item{TN}{The True Negatives (TN).}
#'     \item{FP}{The False Positives (FP).}
#'     \item{FN}{The False Negatives (FN).}
#'     \item{TSS}{The True Skill Score (TSS) calculated as (TP / (TP + FN)) - (FP / (FP + TN)).}
#'     \item{order}{The order of the indicator values after sorting in descending order.}
#'   }
#' 
#' @export
setGeneric("roc", function(state, indicator, ...) {
  standardGeneric("roc")
})

setMethod("roc",
          signature(state = "numeric", indicator = "numeric"),
          function(state, indicator, ...) {
            order = order(indicator, decreasing = TRUE)
            state = state[order]
            indicator = indicator[order]
            label = state > 1
            
            cumn = seq(1, length(label))
            
            TPR = cumsum(label) / sum(label)
            FPR = cumsum(!label) / sum(!label)
            
            TP = cumsum(label)
            FP = cumsum(!label)
            
            TN = sum(!label) - FP
            FN = sum(label) - TP
            
            TSS = (TP / (TP + FN) - FP / (FP + TN))
            
            result_df = data.frame(
              state = state,
              label = label,
              indicator = indicator,
              TPR = TPR,
              FPR = FPR,
              TP = TP,
              TN = TN,
              FP = FP,
              FN = FN,
              TSS = TSS,
              order = order
            )
            
            return(result_df)
          })

#' @examples
#' # In this example, we first generate sample data for state and indicator vectors. 
#' # Generate sample data
#' state <- c(0.5, 2.3, 1.2, 1.8, 3.0, 0.7)
#' indicator <- c(0.6, 2.2, 1.1, 1.9, 2.8, 0.5)
#'
#' # Then, we call the roc function to calculate ROC statistics and print the results.
#' # Calculate ROC statistics
#' roc_result=roc(state, indicator)
#'
#' # Print the ROC statistics
#' roc_result
#'
#' #Plot the ROC curve using the ggplot2 package. 
#' 
#' library(ggplot2)
#' ggplot(roc_result, aes(x = FPR, y = TPR)) +
#'   geom_line() +
#'   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#'   labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)") +
#'   ggtitle("ROC Curve") +
#'   theme_minimal()
