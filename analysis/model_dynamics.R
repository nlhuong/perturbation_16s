library(tidyverse)

# Stick breaking process for Dirichlet Process
stick_break <- function(alpha, n = 100) {
    b <- rbeta(n, 1, alpha)
    p <- numeric(n)
    p[1] <- b[1]
    p[2:n] <- sapply(seq(2, n), function(i) b[i] * prod(1 - b[seq_len(i - 1)]))
    return(p)
}

expand_matrix <- function(c, M) {
    classes <- sort(unique(c))
    if(ncol(M) != nrow(M)) {
        stop("M must be a square matrix.")
    }
    if(ncol(M) < length(classes)) {
        stop("M interaction matrix dimension does not match the number ",
             "of classes")
    }
    n <- length(c)
    edges <- expand.grid(i = seq_len(n), j = seq_len(n)) %>%
        mutate(ci = c[i], cj = c[j])
    mm <- reshape2::melt(M, varnames = c("ci", "cj"))
    df <- suppressMessages(edges %>% left_join(mm))
    M <- df %>% 
        select(i, j, value) %>% 
        spread(j, value) %>%
        select(-i)
}


x_dynamics <- function(x0, t, a1, a2, M, sigma_w) {
    dt_vec <- diff(t)
    M <- as.matrix(M)
    X <- matrix(NA, nrow = length(x0), ncol = length(t))
    X[, 1] <- x0
    for (j in seq(1, length(t) - 1)){
        print(j)
        dt <- dt_vec[j]
        sd <- sqrt(dt) * sigma_w
        mu <- X[, j] + X[, j] * (a0 + a1 * X[, j] + M %*% X[, j]) * dt
        print(summary(mu))
        X[, j+1] = mu
        # X[, j+1] <- sapply(seq_along(mu), function(i) {
        #     rnorm(n = 1, mean = mu[i], sd = sd)})
    }
    return(X)
}

t <- seq(0, 1, length.out = 1000) 
dt <- diff(t)
idx <- c(1, 5, 7, 9, 11, 2, 4, 6, 8, 10, 12, 3, 13)
c <- c(1, 2, 3, 2, 1, 2, 1, 2, 1, 2, 1, 2, 3)
I <- length(c)
B <- matrix(c(0, 3, -1, 0, 0, 0, 2, -4, 0), byrow = TRUE, ncol = 3)
M <- expand_matrix(c, B)
a0 <- rnorm(n = I, sd = 10) 
a1 <- rep(-5, I)
x0 <- floor(runif(5, 15, n = I))
X <- x_dynamics(x0, t, a0, a1, M, sigma_w)


X.long <- reshape2::melt(X, varnames = c("species", "time")) %>%
    filter(species %in% sample(1:I, 10))
ggplot(X.long %>% 
           #filter(species %in% sample(1:I, 10)) %>%
           filter(abs(value) < 100),
       aes(x = time, y = value)) +
    geom_point() + facet_wrap(~species, scales ="free")

# Data -------------------------------------------------------------------------
# i = 1, ..., I: species
# j = 1, ..., J; time
t <- seq(0, 1, by = 0.01)#seq(-20, 20) 
dt <- diff(t)
    
# Hyperparameters --------------------------------------------------------------
I <- 50
J <- length(t)
N <- 5
degfree <- 10
gamma_shape <- 2
p_z <- 0.5

# Higher-level priors ----------------------------------------------------------
alpha <- rgamma(n = 1, shape = gamma_shape)
sigma_a <- 1/rchisq(n = 1, df = degfree)
sigma_b <- 1/rchisq(n = 1, df = degfree)
sigma_w <- 1/rchisq(n = 1, df = degfree)

# For Dirichlet Process --------------------------------------------------------
p_class <- stick_break(alpha)
c <- sample(seq_along(p_class), size = I, 
            prob = p_class, replace = TRUE)
K <- max(c) 
B <- matrix(rnorm(K^2, sd = sigma_b), nrow = K)
diag(B) <- 0
# Generalized Lotka-Voltera dynamics usually have
# sign(bij) = - sign(bji) else species grow indefinitely.
signB <- sign(B)
signB[upper.tri(signB)] <- -t(signB)[upper.tri(signB)]
B <- abs(B) * signB

# Edge selection ---------------------------------------------------------------
Z <- matrix(
    sample(c(0, 1), size = length(B), replace = TRUE, 
    prob = c(p_z, 1 - p_z)),
    nrow = nrow(B))
mm <- B * Z
M <- expand_matrix(c, mm)

# Self Interactions-------------------------------------------------------------
a0 <- rnorm(n = I, sd = sigma_a) 
a1 <- rnorm(n = I, sd = sigma_a)

# Dynamics ---------------------------------------------------------------------
x0 <- rnorm(n = I)

X <- x_dynamics(x0, t, a0, a1, M, sigma_w)

X.long <- reshape2::melt(X, varnames = c("species", "time")) %>%
    filter(species %in% sample(1:I, 10))
ggplot(X.long %>% 
           #filter(species %in% sample(1:I, 10)) %>%
           filter(abs(value) < 100),
       aes(x = time, y = value)) +
    geom_point() + facet_wrap(~species, scales ="free")


# Hyperparameters were set using a technique similar
# to (Bucci et al., 2016), where means of distributions were
# empirically calibrated based on the data and variances were
# set to large values to produce diffuse priors.
# Bucci, Vanni, Tzen, Belinda, Li, Ning, Simmons, Matt,
# Tanoue, Takeshi, Bogart, Elijah, Deng, Luxue, Yeliseyev,
# Vladimir, Delaney, Mary L., Liu, Qing, Olle, Bernat,
# Stein, Richard R., Honda, Kenya, Bry, Lynn, and Gerber,
# Georg K. Mdsine: Microbial dynamical systems
# inference engine for microbiome time-series analyses.
# Genome Biology, 17(1):121, 2016