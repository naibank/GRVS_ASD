library("parallel")
library("data.table")

estimate_power <- function(a, h2_g, h2_gt, n, both_required = FALSE) {
  ##h2_g <- 0.04
  N <- 100000
  ## fuck, there are two parts, argh
  G_mA <- rnorm(N, mean = 0, sd = sqrt(h2_g))
  G_mB <- rnorm(N, mean = 0, sd = sqrt(h2_gt - h2_g))
  G_fA <- rnorm(N, mean = 0, sd = sqrt(h2_g))
  G_fB <- rnorm(N, mean = 0, sd = sqrt(h2_gt - h2_g))
  G_m <- G_mA + G_mB
  G_f <- G_fA + G_fB
  ## noise?
  G_o1A <- rnorm(N, mean = 0, sd = sqrt(h2_g / 2)) + 0.5 * (G_mA + G_fA)
  G_o1B <- rnorm(N, mean = 0, sd = sqrt((h2_gt - h2_g)/ 2)) + 0.5 * (G_mB + G_fB)
  G_o1 <- G_o1A + G_o1B
  ##
  G_o2A <- rnorm(N, mean = 0, sd = sqrt(h2_g / 2)) + 0.5 * (G_mA + G_fA)
  G_o2B <- rnorm(N, mean = 0, sd = sqrt((h2_gt - h2_g)/ 2)) + 0.5 * (G_mB + G_fB)
  G_o2 <- G_o2A + G_o2B
  ## parents
  ##
  eps1 <- rnorm(N, mean = 0, sd = sqrt(1 - h2_gt))
  Y1 <- G_o1 + eps1
  Z1 <- Y1 > qnorm(p = 0.015, mean = 0, sd = 1, lower.tail = FALSE)
  eps2 <- rnorm(N, mean = 0, sd = sqrt(1 - h2_gt))
  Y2 <- G_o2 + eps2
  Z2 <- Y2 > qnorm(p = 0.015, mean = 0, sd = 1, lower.tail = FALSE)
  ## ALSO require unaffected parents?
  if (both_required) {
    who <- which(Z1 & Z2)
  } else {
    who <- which(Z1)
  }
  if (length(who) < n) {
    return(NA)
  }
  who <- who[1:n]
  ## use PRS bits!
  val1 <- (G_o1A[who] - 0.5 * (G_mA[who] + G_fA[who]))
  val2 <- (G_o2A[who] - 0.5 * (G_mA[who] + G_fA[who]))
  o1 <- t.test(val1)
  o2 <- t.test(val2)
  p1 <- o1$p.value
  p2 <- o2$p.value
  m1 <- mean(val1)
  m2 <- mean(val2)
  r <- c(p1 = p1, p2 = p2, m1 = m1, m2 = m2)
  return(r)
  ##return(mean(val))
}

get_power <- function(out, name) {
  out <- out[sapply(out, length) == 4]
  out <- t(sapply(out, I))
  sum(out[, name] <= 0.05, na.rm = TRUE) / sum(is.na(out[, name]) == FALSE)
}

set.seed(134)
#power whole cohort 100%
out1 <- mclapply(1:1000, mc.cores = 4, estimate_power, h2_g = 0.0245, h2_gt = 0.60, n = 255, both_required = FALSE)
get_power(out1, "p1")
get_power(out1, "p2")
#power essential group 91%
out1 <- mclapply(1:1000, mc.cores = 4, estimate_power, h2_g = 0.0245, h2_gt = 0.60, n = 141, both_required = FALSE)
get_power(out1, "p1")

#power equivocal group 47%
out1 <- mclapply(1:1000, mc.cores = 4, estimate_power, h2_g = 0.0245, h2_gt = 0.60, n = 100, both_required = FALSE)
get_power(out1, "p1")

#power complex group 64%
out1 <- mclapply(1:1000, mc.cores = 4, estimate_power, h2_g = 0.0245, h2_gt = 0.60, n = 67, both_required = FALSE)
get_power(out1, "p1")

#power non essential 84%
out1 <- mclapply(1:1000, mc.cores = 4, estimate_power, h2_g = 0.0245, h2_gt = 0.60, n = 114, both_required = FALSE)
get_power(out1, "p1")
