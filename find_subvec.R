#' Find a matching sub-vector
#'
#' Given a vector (`invec`) and a no-larger sub-vector (`subvec`),
#' determine if the latter occurs perfectly.
#' @param invec vector
#' @param subvec vector
#' @return integer positions, length 0 or more
find_subvec <- function(invec, subvec) {
  # Source: https://stackoverflow.com/questions/62377018/how-can-i-find-and-replace-a-specific-sequence-of-numbers-in-a-vector-in-r/62377308#62377308
  sublen <- seq_along(subvec) - 1L
  if (length(subvec) > length(invec)) return(integer(0))
  which(
    sapply(seq_len(length(invec) - length(subvec) + 1L),
           function(i) all(subvec == invec[i + sublen]))
  )
}

test_vec = c(0,0,0,1,1,1,1,0,0,0,1,0,1,0,0,0,0,1,1,1,0,0,1,0,1,0)

sub_vec = c(0,1,0)

repl_vec = c(0,0,0)

indices = find_subvec(test_vec, sub_vec)

for (i in indices) test_vec[i + seq_along(sub_vec) - 1] = repl_vec

