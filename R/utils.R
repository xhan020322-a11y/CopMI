`%||%` <- function(x, y) if (!is.null(x)) x else y

empty_fallback_df <- function() {
  data.frame(method = character(0), stage = character(0), count = integer(0), detail = character(0), stringsAsFactors = FALSE)
}

fallback_record <- function(method, stage, count, detail) {
  if (count == 0L) return(empty_fallback_df())
  data.frame(method = as.character(method), stage = as.character(stage), count = count, detail = as.character(detail), stringsAsFactors = FALSE)
}

combine_fallbacks <- function(...) {
  parts <- Filter(function(x) !is.null(x) && nrow(x) > 0L, list(...))
  if (!length(parts)) return(empty_fallback_df())
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

summarize_messages <- function(x, max_n = 3L) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("")
  shown <- utils::head(x, max_n)
  extra <- length(x) - length(shown)
  paste0(paste(shown, collapse = " | "), if (extra > 0L) paste0(" | ... +", extra, " more") else "")
}
