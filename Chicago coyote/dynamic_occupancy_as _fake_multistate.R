
logit <- function(x) log(x / (1-x))

gamma <- 0.65
eps <- 0.34

gl <- logit(gamma)
el <- logit(eps)

tpm <- diag(4)

# State U
tpm[1,1] <- 3 # to U
tpm[2,1] <- exp(gl) # D
tpm[3,1] <- exp(gl) # N
tpm[4,1] <- exp(gl) # DN

# Check to make sure the probabilities line up
tpm[,1] / sum(tpm[,1])

# State D
tpm[1,2] <- exp(el) * 3 # to U
tpm[2:4,2] <- 1 # all other reference states

# check to make sure the probabilities line up
tpm[,2] / sum(tpm[,2])

# State N
tpm[1,3] <- exp(el) * 3 # to U
tpm[2:4,3] <- 1 # all other reference states

# check to make sure the probabilities line up
tpm[,3] / sum(tpm[,3])

# State DN
tpm[1,4] <- exp(el) * 3
tpm[2:4,4] <- 1 # all other reference states

tpm[,4] / sum(tpm[,4])
