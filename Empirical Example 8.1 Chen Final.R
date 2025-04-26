##################################################
# 0. Απαραίτητες Ρυθμίσεις & Βοηθητικές Συναρτήσεις
##################################################
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off(); rm(list = ls(all = TRUE)); cat("\014"); getwd()
# ----------------------------------------
# Συνάρτηση logistic_function (αντιστοιχεί στο logistic.m)
# ----------------------------------------


# ----------------------------------------
# Συνάρτηση BaiNg_selectFactors (ανάλογο του BaiNg.m)
# ----------------------------------------
BaiNg_selectFactors <- function(X, Rmax , c){
  # X: matrix T x N
  Tn <- nrow(X)
  Nn <- ncol(X)
  Smax <- Rmax
  
  # υπολογίζουμε τον πίνακα συνδιασποράς X'X/(N*T)
  covX <- (1/(Nn*Tn)) * crossprod(X)  # (N x N)
  e_decomp <- eigen(covX, symmetric = TRUE)
  eigvals   <- Re(e_decomp$values)
  idx_sort  <- order(eigvals, decreasing=TRUE)
  lambda_sorted <- eigvals[idx_sort]
  
  # Yπολογίζουμε το κριτήριο IC (1,2,3) για k=1..Rmax
  LIC1 <- numeric(Smax)
  LIC2 <- numeric(Smax)
  LIC3 <- numeric(Smax)
  
  sum_all <- sum(lambda_sorted)   # άθροισμα όλων των ιδιοτιμών
  for(k in 1:Smax){
    # reconstruction error ~ sum_all - sum(lambda_sorted[1:k])
    re <- sum_all - sum(lambda_sorted[1:k])
    # Τύποι IC κατά Bai & Ng (2002), π.χ.:
    LIC1[k] <- log(re) + c * k * (Nn + Tn)/(Nn*Tn)* log((Nn*Tn)/(Nn+Tn))
    LIC2[k] <- log(re) + c * k * (Nn + Tn)/(Nn*Tn)* log(min(Nn, Tn))
    LIC3[k] <- log(re) + c * k * log(min(Nn, Tn))/min(Nn, Tn)
  }
  IC1 <- which.min(LIC1)
  IC2 <- which.min(LIC2)
  IC3 <- which.min(LIC3)
  
  return(c(IC1, IC2, IC3))
}


##################################################
# 1. Φόρτωση & Κανονικοποίηση Δεδομένων
##################################################
library(readr)
library(GrFA)

# Προσαρμόστε τα paths στα δικά σας αρχεία:
file_dev  <- "dev1.csv"
file_eme  <- "eme1.csv"

library(readxl)
data_dev  <- read_csv(file_dev)
data_eme  <- read_csv(file_eme)
# π.χ. Q3-98

# Μετατρέπουμε σε απλούς πίνακες
X_dev <- as.matrix(data_dev)
X_eme <- as.matrix(data_eme)

Tn <- nrow(X_dev)   # Υποθέτουμε ίδιο πλήθος γραμμών = χρόνοι


# Ενώνουμε
X0 <- cbind(X_dev, X_eme)  # (Tn x (N1+N2))

# Standardization (στήλη-στήλη)
Xmean <- colMeans(X0)
Xsd   <- apply(X0, 2, sd)
X     <- scale(X0, center = Xmean, scale = Xsd)  # Tn x (N1+N2)

N1 <- ncol(X_dev)
N2 <- ncol(X_eme)
##################################################
# 2. Εκτίμηση Τοπικών (Within-Group) Παραγόντων
##################################################
# Έστω M=2 (Group1=Developed, Group2=Emerging)
M <- 2
Smax <- 15  # το πολύ 15 factors ανά ομάδα, όπως στο Matlab example

X1 <- X[, 1:N1, drop=FALSE]           # Tn x N1
X2 <- X[, (N1+1):(N1+N2), drop=FALSE] # Tn x N2
X3 <- as.matrix(X1)
X4 <- as.matrix(X2)
anyNA(X3)
anyNA(X4)
i_values <- 1:floor(N1/2)

# Αρχικοποίηση κενών vectors για τις τιμές IC2 και ER
est_IC2 <- c()
est_ER <- c()

# Υπολογισμός τιμών για τις στήλες est_num(...)
for (i in i_values) {
  # Υπολογισμός των τιμών
  # IC2_value <- suppressWarnings(est_num(X3, kmax = i, type = "IC2"))
  # ER_value <- suppressWarnings(est_num(X3, kmax = i, type = "ER"))
  IC2_value <- est_num(X3, kmax = i, type = "IC2")
  ER_value <- est_num(X3, kmax = i, type = "ER")
  # Προσθήκη των τιμών μόνο αν δεν είναι NA
  if (!is.na(IC2_value)) {
    est_IC2 <- c(est_IC2, IC2_value)
  }
  if (!is.na(ER_value)) {
    est_ER <- c(est_ER, ER_value)
  }
}

# Δημιουργία του data.frame
results1 <- data.frame(
  i = i_values,
  est_IC2,
  est_ER
)

# Εμφάνιση του αποτελέσματος
print(results1)
i_values <- 1:floor(N2/2)

# Αρχικοποίηση κενών vectors για τις τιμές IC2 και ER
est_IC2 <- c()
est_ER <- c()

# Υπολογισμός τιμών για τις στήλες est_num(...)
for (i in i_values) {
  # Υπολογισμός των τιμών
  IC2_value <- est_num(X4, kmax = i, type = "IC2")
  ER_value <- est_num(X4, kmax = i, type = "ER")
  
  # Προσθήκη των τιμών μόνο αν δεν είναι NA
  if (!is.na(IC2_value)) {
    est_IC2 <- c(est_IC2, IC2_value)
  }
  if (!is.na(ER_value)) {
    est_ER <- c(est_ER, ER_value)
  }
}

# Δημιουργία του data.frame
results2 <- data.frame(
  i = i_values,
  est_IC2,
  est_ER
)

# Εμφάνιση του αποτελέσματος
print(results2)
# Ορισμός των labels
# Ορισμός των labels για τις εκτιμήσεις
# Ορισμός των labels για τις εκτιμήσεις

array <- c("PC1", "PC2", "PC3", "IC1", "IC2", "IC3", "AIC3", "BIC3", "ER", "GR")

# Ορισμός των εύρων τιμών i
i_values <- 1:floor(N1/2)

# Δημιουργία κενού πίνακα δεδομένων με 36 γραμμές και 11 στήλες
results1 <- data.frame(matrix(ncol = length(array) + 1, nrow = length(i_values)))

# Ορισμός ονομάτων στηλών (συμπεριλαμβανομένου του i)
colnames(results1) <- c("i", paste0("est_", array))

# Συμπλήρωση της στήλης i με τις αντίστοιχες τιμές
results1$i <- i_values

# Υπολογισμός των τιμών για κάθε i και αποθήκευση στο data.frame
for (i in i_values) {
  for (j in seq_along(array)) {
    # Υπολογισμός της εκτίμησης
    # value <- suppressWarnings(est_num(X3, kmax = i, type = array[j]))
    value <- est_num(X3, kmax = i, type = array[j])
    
    # Αν δεν είναι NA, αποθήκευση της τιμής
    results1[i, j + 1] <- ifelse(!is.na(value), value, NA)
  }
}

# Εκτύπωση του αποτελέσματος
print(results1)

array <- c("PC1", "PC2", "PC3", "IC1", "IC2", "IC3", "AIC3", "BIC3", "ER", "GR")

# Ορισμός των εύρων τιμών i
i_values <- 1:floor(N2/2)

# Δημιουργία κενού πίνακα δεδομένων με 36 γραμμές και 11 στήλες
results2 <- data.frame(matrix(ncol = length(array) + 1, nrow = length(i_values)))

# Ορισμός ονομάτων στηλών (συμπεριλαμβανομένου του i)
colnames(results2) <- c("i", paste0("est_", array))

# Συμπλήρωση της στήλης i με τις αντίστοιχες τιμές
results2$i <- i_values

# Υπολογισμός των τιμών για κάθε i και αποθήκευση στο data.frame
for (i in i_values) {
  for (j in seq_along(array)) {
    # Υπολογισμός της εκτίμησης
    value <- est_num(X4, kmax = i, type = array[j])
    
    # Αν δεν είναι NA, αποθήκευση της τιμής
    results2[i, j + 1] <- ifelse(!is.na(value), value, NA)
  }
}

# Εκτύπωση του αποτελέσματος
print(results2)


# Εκτύπωση του αποτελέσματος

kkk <- 0
M <- 2
T <- nrow(X)
N <- c(N1, N2)
Proj <- array(0, dim = c(T, T, M))
CirProjection <- diag(T)
WithinGroupFactorsRatioTest <- rep(0, M)
WithinGroupFactorsIC <- array(0, dim = c(1, 3, M))
IC2 <- rep(0, M)

for (s in 1:M) {
  Nm <- N[s]
  Xm <- X[, (kkk+1):(kkk+Nm)]
  eig_result <- eigen(1/(Nm*T) * t(Xm) %*% Xm)
  vector <- eig_result$vectors
  value <- eig_result$values
  
  # Έλεγχος αν υπάρχουν επαρκείς ιδιοτιμές πριν την ταξινόμηση
  if (length(value) < 10) {
    stop("Error: Not enough eigenvalues for sorting.")
  }
  
  sort_result <- sort(Re(value), decreasing = TRUE)
  
  # Έλεγχος μήκους του sort_result για αποφυγή out-of-bounds access
  if (length(sort_result) < 10) {
    stop("Error: sort_result has fewer than 10 elements.")
  }
  
  position <- which.max(sort_result[1:9] / sort_result[2:10])
  WithinGroupFactorsRatioTest[s] <- position
  
  Rmax <- min(15, length(sort_result))  # Περιορισμός στο διαθέσιμο μέγεθος
  c <- 1
  WithinGroupFactorsIC[1, 1:3, s] <- BaiNg_selectFactors(Xm, Rmax, c)
  IC2[s] <- WithinGroupFactorsIC[1, 2, s]
  
  Llambda <- matrix(0, nrow = Nm, ncol = Smax)
  Lfhat <- matrix(0, nrow = T, ncol = Smax)
  LChat <- matrix(0, nrow = T, ncol = Nm)
  
  Sm <- min(15, ncol(vector))  # Προσαρμογή ώστε να μην υπερβαίνει το μήκος των δεδομένων
  
  if (Sm > 0) {
    Llambda[1:Nm, 1:Sm] <- sqrt(Nm) * Re(vector[, 1:Sm])
    Lfhat[1:T, 1:Sm] <- (1/Nm) * Xm %*% Llambda[1:Nm, 1:Sm]
    LChat[1:T, 1:Nm] <- Lfhat[1:T, 1:Sm] %*% t(Llambda[1:Nm, 1:Sm])
  } else {
    warning("Sm is zero, skipping Llambda, Lfhat, and LChat calculations.")
  }
  
  Ldevi <- norm(LChat[1:T, 1:Nm] - Xm, "F") / sqrt(Nm * T)
  
  if (Sm > 0) {
    KmHat <- Lfhat[1:T, 1:Sm]
  } else {
    KmHat <- matrix(0, nrow = T, ncol = 1)
  }
  
  if (s == 1) {
    K1hat <- KmHat
  }
  
  if (ncol(KmHat) > 0 && nrow(KmHat) > 0) {
    Proj[1:T, 1:T, s] <- KmHat %*% MASS::ginv(t(KmHat) %*% KmHat) %*% t(KmHat)
  } else {
    Proj[1:T, 1:T, s] <- diag(T)
  }
  
  CirProjection <- CirProjection %*% Proj[1:T, 1:T, s]
  
  kkk <- kkk + N[s]
}
get_pca_factors <- function(Xg, s=15){
  Tn <- nrow(Xg)
  Ng <- ncol(Xg)
  covXg <- (1/(Ng*Tn)) * crossprod(Xg)
  eig_out <- eigen(covXg, symmetric=TRUE)
  vals  <- Re(eig_out$values)
  vecs  <- Re(eig_out$vectors)
  idx_sort <- order(vals, decreasing=TRUE)
  vals_sorted <- vals[idx_sort]
  vecs_sorted <- vecs[, idx_sort, drop=FALSE]
  
  s_use <- min(s, length(vals_sorted))
  Lambda <- sqrt(Ng) * vecs_sorted[, 1:s_use]
  # Khat = (1/Ng) * Xg %*% Lambda
  Khat <- (1/Ng) * Xg %*% Lambda   # (Tn x s_use)
  return(Khat)
}

Khat1 <- get_pca_factors(X1, s=Smax)   # Tn x s1
Khat2 <- get_pca_factors(X2, s=Smax)   # Tn x s2
s1 <- ncol(Khat1)
s2 <- ncol(Khat2)


##################################################
# 3. Circular Projection για Global Factors
##################################################
# Πίνακες προβολής P(Khat_m):
projK1 <- Khat1 %*% solve(crossprod(Khat1)) %*% t(Khat1)  # (Tn x Tn)
projK2 <- Khat2 %*% solve(crossprod(Khat2)) %*% t(Khat2)

# Προϊόν π.χ. P(K2)*P(K1) (για M=2)
tmpProd <- projK2 %*% projK1
Upsilon_hat <- t(tmpProd) %*% tmpProd  # (Tn x Tn)

# Ιδιοδιάσπαση της Upsilon_hat
ev_U <- eigen(Upsilon_hat, symmetric=TRUE)
valsU <- Re(ev_U$values)
vecsU <- Re(ev_U$vectors)
idxU  <- order(valsU, decreasing=TRUE)
valsU_sorted <- valsU[idxU]
vecsU_sorted <- vecsU[, idxU]

# Πόσους global factors θες; π.χ. R=5
R <- 5  
globalEvecs <- vecsU_sorted[, 1:R, drop=FALSE]  # (Tn x R)

# Εκτιμώμενοι global factors
Fhat_circ <- sqrt(Tn) * globalEvecs   # (Tn x R)


##################################################
# 4. Augmented Circular Projection
#    (reference group = Developed = group1)
##################################################
Km0 <- Khat1
# (Khat1'Khat1)^{-1/2}
# Μια γρήγορη λύση είναι chol->solve, ή μπορείς να χρησιμοποιήσεις svd/eigen.
invKm0 <- solve(chol(crossprod(Km0)))

tmpProdAug   <- projK2 %*% projK1  # (Tn x Tn)
middle_part  <- t(tmpProdAug) %*% tmpProdAug   # (Tn x Tn)

left_part  <- invKm0 %*% t(Km0)  # (s1 x Tn)
right_part <- Km0 %*% invKm0     # (Tn x s1)

Ahat <- left_part %*% middle_part %*% right_part  # (s1 x s1)

evA <- eigen(Ahat, symmetric=TRUE)
valsA <- Re(evA$values)
vecsA <- Re(evA$vectors)
idxA  <- order(valsA, decreasing=TRUE)
valsA_sorted <- valsA[idxA]
vecsA_sorted <- vecsA[, idxA, drop=FALSE]
Raug <- 5
phi_hat <- vecsA_sorted[, 1:Raug, drop=FALSE]  # (s1 x Raug)
Fhat_aug <- sqrt(Tn) * Km0 %*% invKm0 %*% phi_hat  # (Tn x Raug)

##################################################
# 5. Παράδειγμα Table_4 & Table_4(aug)
##################################################
# Π.χ. υπολογίζουμε διακύμανση του πρώτου τοπικού παράγοντα vs πρώτου global.
Lfhat1_dev <- Khat1[, 1, drop=FALSE]
Lfhat1_eme <- Khat2[, 1, drop=FALSE]
GF1   <- Fhat_circ[,1,drop=FALSE]
AGF1  <- Fhat_aug[,1,drop = FALSE]
s <- 1

# Υπολογισμός SVD για K1hat' * K1hat
svd_result <- svd(t(K1hat) %*% K1hat)
Uhat <- svd_result$u
Dhat <- diag(svd_result$d)  # Μετατροπή των singular values σε διαγώνιο πίνακα
Vhat <- svd_result$v

# Υπολογισμός της ημι-δύναμης του K1hat' * K1hat
K1K1hathalf <- Uhat %*% sqrt(Dhat) %*% t(Uhat)

# Υπολογισμός SVD για προβολή
svd_hat <- svd((t(CirProjection %*% K1hat %*% solve(K1K1hathalf))) %*% CirProjection %*% K1hat %*% solve(K1K1hathalf))
Usvdhat <- svd_hat$u
Dsvdhat <- diag(svd_hat$d)
Vsvdhat <- svd_hat$v

# Ορισμός s_star
s_star <- Sm

# Υπολογισμός προβολικού σκορ
regresScore <- sort(
  diag(t(K1hat %*% solve(K1K1hathalf) %*% Usvdhat[, 1:s_star]) %*% 
         (M * diag(T) - apply(Proj[, , 1:M], c(1, 2), sum)) %*% 
         (K1hat %*% solve(K1K1hathalf) %*% Usvdhat[, 1:s_star])) / M, 
  decreasing = FALSE
)

# Ορισμός λογιστικής συνάρτησης
logistic <- function(P1, P0, tau, x) {
  A <- P1 / P0 - 1
  P1 / (1 + A * exp(-tau * x))
}

# Υπολογισμός αναλογίας προβολής
projectionRatio <- logistic(1, 0.001, 14, log(log(min(T, Nm))) * regresScore[2:s_star]) - 
  logistic(1, 0.001, 14, log(log(min(T, Nm))) * regresScore[1:(s_star-1)])

# Εύρεση μέγιστης τιμής και θέσης
maxValue1 <- max(projectionRatio)
maxPosition1 <- which.max(projectionRatio)
globalNoRatio <- maxPosition1

# Καθορισμός αριθμού global factors
R <- globalNoRatio
phi <- Usvdhat[, 1:R]

# Έλεγχος αν το phi είναι ορθοκανονικό
t(phi) %*% phi

# Υπολογισμός παγκόσμιων παραγόντων
Fhat <- sqrt(T) * (K1hat %*% solve(K1K1hathalf) %*% phi)

# Υπολογισμός φορτίων
Lambdahat <- 1/T * t(X) %*% Fhat

# Υπολογισμός κοινών global components
globalCommon <- Fhat %*% t(Lambdahat)
ACP <- function(y, reference_group, rmax = 8, r0 = NULL, r = NULL, r_aug = NULL, localfactor = FALSE, type = "IC3") {
  M <- length(y)
  T <- nrow(y[[1]])
  Nm <- sapply(y, ncol)
  K <- list()
  Proj_Mat <- list()
  
  for (m in 1:M) {
    K[[m]] <- FA(y[[m]], r = rmax)$F  # T x rmax
    Proj_Mat[[m]] <- K[[m]] %*% solve(t(K[[m]]) %*% K[[m]]) %*% t(K[[m]])  # T x T
  }
  
  Aug_Proj_Mat <- Proj_Mat[[reference_group]]
  for (m in setdiff(1:M, reference_group)) {
    Aug_Proj_Mat <- Aug_Proj_Mat %*% Proj_Mat[[m]]
  }
  
  Aug_Proj_Mat <- t(Aug_Proj_Mat) %*% Aug_Proj_Mat
  eig_deco <- eigen(Aug_Proj_Mat)
  eig_vec <- eig_deco$vectors[, 1:rmax, drop = FALSE]
  rho <- eig_deco$values[1:rmax]
  
  if (is.null(r0)) {
    ARSS <- rep(0, rmax)
    for (i in 1:rmax) {
      for (m in 1:M) {
        ARSS[i] <- ARSS[i] + t(eig_vec[, i]) %*% (diag(T) - Proj_Mat[[m]]) %*% eig_vec[, i]
      }
    }
    ARSS <- ARSS / M
    nmin <- min(c(T, Nm))
    logistic <- function(x) {
      P1 <- 1
      P0 <- 10^-3
      A <- P1 / P0 - 1
      tau <- 14
      P1 / (1 + A * exp(-tau * x))
    }
    r0hat <- which.max(logistic(log(log(nmin)) * ARSS[2:rmax]) - logistic(log(log(nmin)) * ARSS[1:(rmax - 1)]))
  } else {
    if (!(r0 %% 1 == 0) | r0 <= 0) stop("invalid 'r0' input")
    r0hat <- r0
  }
  
  Ghat <- eig_deco$vectors[, 1:r0hat, drop = FALSE] * sqrt(T)  # T x r0hat
  Proj_G <- Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat)        # T x T
  
  loading_G <- list()
  for (m in 1:M) {
    loading_G[[m]] <- t(y[[m]]) %*% Ghat / T  # Nm x r0hat
  }
  
  eig_aug <- eigen(Aug_Proj_Mat)
  if (is.null(r_aug)) {
    r_aug <- sum(eig_aug$values > 1e-4)
  } else {
    if (!(r_aug %% 1 == 0) | r_aug <= 0) stop("invalid 'r_aug' input")
  }
  
  AGhat <- eig_aug$vectors[, 1:r_aug, drop = FALSE] * sqrt(T)  # T x r_aug
  loading_AG <- list()
  for (m in 1:M) {
    loading_AG[[m]] <- t(y[[m]]) %*% AGhat / T  # Nm x r_aug
  }
  
  rhat <- rep(0, M)
  
  if (localfactor == FALSE) {
    res <- list(
      r0hat = r0hat,
      rhat = rhat,
      r_aug = r_aug,
      rho = rho,
      Ghat = Ghat,
      AGhat = AGhat,
      loading_G = loading_G,
      loading_AG = loading_AG
    )
  } else {
    Fhat <- list()
    loading_F <- list()
    y_proj_G <- lapply(y, function(x) x - Proj_G %*% x)
    
    if (is.null(r)) {
      for (m in 1:M) {
        rhat[m] <- est_num(y_proj_G[[m]], kmax = rmax - r0hat, type = type)
        fit <- FA(y_proj_G[[m]], r = rhat[m])
        Fhat[[m]] <- fit$F  # T x rhat
        loading_F[[m]] <- fit$L  # Nm x rhat
      }
    } else {
      if (!(all(r %% 1 == 0) && all(r >= 0))) stop("invalid r input")
      rhat <- r
      for (m in 1:M) {
        fit <- FA(y_proj_G[[m]], r = rhat[m])
        Fhat[[m]] <- fit$F
        loading_F[[m]] <- fit$L
      }
    }
    
    e <- list()
    for (m in 1:M) {
      if (rhat[m] > 0) {
        e[[m]] <- y_proj_G[[m]] - Fhat[[m]] %*% t(loading_F[[m]])
      } else {
        e[[m]] <- y_proj_G[[m]]
      }
    }
    
    res <- list(
      r0hat = r0hat,
      rhat = rhat,
      r_aug = r_aug,
      rho = rho,
      Ghat = Ghat,
      AGhat = AGhat,
      Fhat = Fhat,
      loading_G = loading_G,
      loading_AG = loading_AG,
      loading_F = loading_F,
      residual = e
    )
  }
  
  cat("The number of global factors is", r0hat, "\n")
  cat("The number of local factors is", rhat, "\n")
  cat("The number of augmented global factors is", r_aug, "\n")
  
  class(res) <- "GFA"
  return(res)
}

out3 <- ACP(list(X1,X2),reference_group = 1,rmax = 15,r0 = NULL,r = NULL,r_aug = NULL,localfactor = TRUE,type = "ER")
# Αποθήκευση αποτελεσμάτων σε CSV αρχεία
write.csv(globalCommon, "globalCommonComponentsAug.csv", row.names = FALSE)
write.csv(out3$AGhat, "globalFactorsAug.csv", row.names = FALSE)
write.csv(Fhat,"globalFactors.csv",row.names = FALSE)
kkk <- 0
LChat2 <- array(0, dim = c(T, sum(N)))

WithinVariations <- numeric(M)
GlobalVariations <- numeric(M)
devi <- numeric(M)


for (s in 1:M) {
  Nm <- N[s]
  Xm <- X[, (kkk + 1):(kkk + Nm)]
  
  eig_result <- eigen(1/(Nm*T) * t(Xm) %*% Xm)
  sort_indices <- order(Re(eig_result$values), decreasing = TRUE)
  Sm <- 1
  
  # Εξαγωγή ιδιοδιανυσμάτων
  loading_s <- sqrt(Nm) * Re(eig_result$vectors[, sort_indices[1:Sm]])
  Fhat_s <- 1/Nm * Xm %*% loading_s
  
  out3$loading_F[[s]] <- loading_s
  out3$Fhat[[s]] <- Fhat_s
  
  LChat2[1:T, (kkk + 1):(kkk + Nm)] <- Fhat_s %*% t(loading_s)
  
  devi[s] <- norm(LChat2[1:T, (kkk + 1):(kkk + Nm)] - Xm, 'f') / sqrt(Nm*T)
  WithinVariations[s] <- mean(apply(LChat2[1:T, (kkk + 1):(kkk + Nm)], 2, var))
  
  # Έλεγχος για AGhat και loading_AG
  
  if (exists("globalCommon")) {
    GlobalVariations[s] <- mean(apply(globalCommon[1:T, (kkk + 1):(kkk + Nm)], 2, var))
  } else {
    GlobalVariations[s] <- NA
  }
  
  
  
  kkk <- kkk + Nm
}

Table4_aug <- round(cbind(WithinVariations, GlobalVariations), 3)
print(Table4_aug)


write.csv(cbind(WithinVariations, GlobalVariations), 'Table 4(aug).csv', row.names = FALSE)
out <- CP(list(X1, X2), rmax = 15, r0 = NULL,r = NULL, localfactor = TRUE, type = "ER")
Lambdahat <- as.matrix(out$loading_G[[1]])
globalCommon <- Fhat %*% t(Lambdahat)
Lambdahat <- as.matrix(out$loading_G[[2]])
globalCommon <- Fhat %*% t(Lambdahat)
Lambdahat <- as.matrix(do.call(rbind, out$loading_G))
globalCommon <- Fhat %*% t(Lambdahat)

M <- 2
kkk <- 0
devi <- numeric(M)
WithinVariations <- numeric(M)
GlobalVariations <- numeric(M)
LChat2 <- matrix(0, nrow = T, ncol = ncol(X))

for (s in 1:M) {
  Nm <- N[s]
  Xm <- X[, (kkk + 1):(kkk + Nm)]
  
  # Υπολογισμός τοπικού παραγοντικού μοντέλου (Sm=1)
  eigen_result <- eigen((1 / (Nm * T)) * t(Xm) %*% Xm)
  sort_result <- order(Re(eigen_result$values), decreasing = TRUE)
  Sm <- 1
  
  out$loading_F[[s]] <- sqrt(Nm) * Re(eigen_result$vectors[, sort_result[1:Sm]])
  out$Fhat[[s]] <- 1 / Nm * Xm %*% out$loading_F[[s]]
  LChat2[, (kkk + 1):(kkk + Nm)] <- out$Fhat[[s]] %*% t(out$loading_F[[s]])
  
  WithinVariations[s] = mean(apply(LChat2[1:T, (kkk + 1):(kkk + Nm)], 2, var))
  
  if (exists("globalCommon")) {  # Έλεγχος αν υπάρχει το globalCommon
    GlobalVariations[s] = mean(apply(globalCommon[1:T, (kkk + 1):(kkk + Nm)], 2, var))
  } else {
    GlobalVariations[s] = NA  # Αν δεν υπάρχει, αποθηκεύουμε NA
  }
  
  kkk = kkk + Nm
}

Table4 <- round(cbind(WithinVariations, GlobalVariations), 3)
print(Table4)
write.csv(Table4, "Table4.csv", row.names = FALSE)
write.csv(cbind(WithinVariations, GlobalVariations), 'Table 4.csv', row.names = FALSE)

# Φόρτωση των δεδομένων από τα αρχεία CSV
table1 <- read.csv("Table 4.csv")
table2 <- read.csv("Table 4(aug).csv")

# Δημιουργία του συνδυασμένου πίνακα με τις επιθυμητές τιμές
combined_table <- rbind(
  cbind(table1[1, 1], table1[1, 2], table2[1, 2]),
  cbind(table1[2, 1], table1[2, 2], table2[2, 2])
)

# Ορισμός των ονομάτων των στηλών
colnames(combined_table) <- c("The first within-group", "Global", "Aug.Global")

# Ορισμός των ονομάτων των γραμμών
rownames(combined_table) <- c("Developed (group 1):", "Emerging (group 2):")

# Στρογγυλοποίηση των τιμών στο 3ο δεκαδικό ψηφίο
combined_table <- round(combined_table, 3)

# Εκτύπωση του τελικού πίνακα
print(combined_table)



library(readr)
library(scales)
X <- (X0 - matrix(rep(colMeans(X0), each = nrow(X0)), nrow = nrow(X0), byrow = FALSE)) / matrix(rep(apply(X0, 2, sd), each = nrow(X0)), nrow = nrow(X0), byrow = FALSE)
Fhataug <- read_csv("globalFactorsAug.csv")  # Augmented global factors
Fhat <- read_csv("globalFactors.csv")        # Global factors



# Επιλογή δεδομένων για τα γραφήματα
X_points1 <- X1[, 10]  
X_points2 <- X2[, 50]
X_points3 <- X[, 50]
X_points4 <- X[, 10]
# Ορισμός ετικετών και δεικτών χρόνου
time_labels <- c("Q3/98", "Q2/01", "Q1/04", "Q4/06", "Q3/09", "Q2/12", "Q1/15", "Q4/17")
time_indices <- seq(1, length(Fhat[[1]]), length.out = length(time_labels))

# Προσαρμογή του υπομνήματος
add_legend <- function(...) {
  legend("topright", inset = 0.03, cex = 0.6, ...)
}

# Συνάρτηση για τη δημιουργία γραφημάτων (με προαιρετικό υπόμνημα)
plot_graph <- function(y1, y2, points, title, y_ticks, legend_labels = NULL) {
  plot(y1, type = "l", col = "blue", lwd = 2, ylim = range(y_ticks), xaxt = "n", yaxt = "n",
       xlab = "Χρόνος", ylab = "Τιμή", main = title)
  lines(y2, col = "red", lwd = 2, lty = 2)
  points(1:length(points), points, col = "black", pch = 21, bg = "black", cex = 1)
  axis(1, at = time_indices, labels = time_labels, cex.axis = 0.8)
  axis(2, at = y_ticks, labels = y_ticks, las = 1, cex.axis = 0.8)
  
  # Προσθήκη υπομνήματος μόνο αν έχει δοθεί
  if (!is.null(legend_labels)) {
    add_legend(legend = legend_labels, col = c("blue", "red"), lty = c(1,2), lwd = 2)
  }
}

# Ορισμός τιμών για άξονα y
y_ticks1 <- seq(-6, 4, by = 2)
y_ticks2 <- seq(-6, 2, by = 1)

# Δημιουργία διαγραμμάτων
dev.new()
par(mfrow = c(3, 2))  # 3 γραμμές, 2 στήλες

# (a) Global vs. Within (Emerging) (με υπόμνημα)
plot_graph(Fhat[[1]], -Lfhat1_eme, X_points2, "(a) Global vs. Within (Emerging)", y_ticks1, 
           c("global factor", "within factor (emerging)"))

# (b) Global vs. Within (Developed) (Χωρίς υπόμνημα)
plot_graph(Fhat[[1]], -Lfhat1_dev, X_points1, "(b) Global vs. Within (Developed)", y_ticks2)

# (c) Aug. Global vs. Within (Emerging) (με υπόμνημα)
plot_graph(-Fhataug[[1]], -Lfhat1_eme, X_points2, "(c) Aug. Global vs. Within (Emerging)", y_ticks1, 
           c("aug. global factor", "within factor (emerging)"))

# (d) Aug. Global vs. Within (Developed) (Χωρίς υπόμνημα)
plot_graph(-Fhataug[[1]], -Lfhat1_dev, X_points1, "(d) Aug. Global vs. Within (Developed)", y_ticks2)

# (e) Global vs. Aug. Global (Emerging) (με υπόμνημα)
plot_graph(Fhat[[1]], -Fhataug[[1]], X_points2, "(e) Global vs. Aug. Global (Emerging)", y_ticks1, 
           c("global factor", "aug. global factor"))

# (f) Global vs. Aug. Global (Developed) (Χωρίς υπόμνημα)
plot_graph(Fhat[[1]], -Fhataug[[1]], X_points1, "(f) Global vs. Aug. Global (Developed)", y_ticks2)
dev.new()
par(mfrow = c(3,2))# (a) Global vs. Within (Emerging) (με υπόμνημα)
plot_graph(Fhat[[1]], -Lfhat1_eme, X_points3, "(a) Global vs. Within (Emerging)", y_ticks1, 
           c("global factor", "within factor (emerging)"))

# (b) Global vs. Within (Developed) (Χωρίς υπόμνημα)
plot_graph(Fhat[[1]], -Lfhat1_dev, X_points4, "(b) Global vs. Within (Developed)", y_ticks2)

# (c) Aug. Global vs. Within (Emerging) (με υπόμνημα)
plot_graph(-Fhataug[[1]], -Lfhat1_eme, X_points3, "(c) Aug. Global vs. Within (Emerging)", y_ticks1, 
           c("aug. global factor", "within factor (emerging)"))

# (d) Aug. Global vs. Within (Developed) (Χωρίς υπόμνημα)
plot_graph(-Fhataug[[1]], -Lfhat1_dev, X_points4, "(d) Aug. Global vs. Within (Developed)", y_ticks2)

# (e) Global vs. Aug. Global (Emerging) (με υπόμνημα)
plot_graph(Fhat[[1]], -Fhataug[[1]], X_points3, "(e) Global vs. Aug. Global (Emerging)", y_ticks1, 
           c("global factor", "aug. global factor"))

# (f) Global vs. Aug. Global (Developed) (Χωρίς υπόμνημα)
plot_graph(Fhat[[1]], -Fhataug[[1]], X_points4, "(f) Global vs. Aug. Global (Developed)", y_ticks2)


