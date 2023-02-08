vec_size <- 5L
n_samp <- 2L
my_samps <- replicate(n = 20000, expr = {
  idxs <- seq(1L, vec_size)
  s1 <- sample(x = idxs, size = n_samp, replace = FALSE)
  remaining <- idxs[-s1]
  
  new_idx_sample_prob <- (n_samp + 1)/(vec_size + 1)
  length_remaining <- vec_size - n_samp
  
  s2 <- sample(x = c(remaining, vec_size + 1L), size = 1,
               prob = c(rep((1-new_idx_sample_prob)/length_remaining, length_remaining), new_idx_sample_prob))
  c(s1, s2)
}, simplify = FALSE)

sapply(my_samps, function(samp) {
  6L %in% samp[1]
}) |> mean()


##################################
##################################

n_treatment <- 10L
n_control <- 100L
n_total <- n_treatment + n_control

t = 0
m = 0
samples <- integer(length = n_treatment)

while (m < n_treatment) {
  u = runif(1)
  if ((n_total - t) * u >= n_treatment - m) {
    t <- t + 1
  } else {
    samples[m + 1] <- t
    t <- t + 1
    m <- m + 1
  }
}

# a strategy for sampling random numbers, over the range (0, ..., k-1, k), where k gets mass p, and mass 1-p is 

k <- 6
p <- 0.8
u <- runif(200000)
l <- u > (1 - p)
ifelse(l, k, floor(u/(1 - p) * k))


###############
N_c <- 10
N_t <- 5

samps <- replicate(n = 10000, expr = {
  a <- seq(1L, N_c + 1)
  v <- integer(length = N_t)
  for (k in 1:N_t) {
    p <- k / (N_c + k)
    u <- runif(1)
    # obtain the position to sample
    # u <- runif(100000)
    samp_pos <- ifelse(u > (1 - p), N_c + 1, floor(u/(1 - p) * N_c) + 1)
    # sample at that position
    v[k] <- a[samp_pos]
    # swap the final position with samp_pos
    a[samp_pos] <- a[N_c + 1]
    # put next treatment cell into rightmost idx of a
    a[N_c + 1] <- k + N_c + 1
  }
  v
}, simplify = FALSE)

# verify key inductive WOR sampling property empirically
sapply(X = samps, FUN = function(samp) {
  any(12 == samp[1:2])
}) |> mean()
