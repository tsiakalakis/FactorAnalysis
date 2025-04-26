rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off(); rm(list = ls(all = TRUE)); cat("\014"); getwd()
# ----------------------------------------
# Συνάρτηση logistic_function (αντιστοιχεί στο logistic.m)
# ----------------------------------------


# ----------------------------------------
# Συνάρτηση BaiNg_selectFactors (ανάλογο του BaiNg.m)
# ----------------------------------------



##################################################
# 1. Φόρτωση & Κανονικοποίηση Δεδομένων
##################################################
library(readr)
library(GrFA)

# Προσαρμόστε τα paths στα δικά σας αρχεία:
file_dev  <- "dev1.csv"
file_eme  <- "eme1.csv"

library(readxl)
data_dev  <- read.csv(file_dev,header = F)
data_eme  <- read.csv(file_eme,header = F)
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
array2 <- c("IC1","IC2","IC3")
  for (s in 1:M){
    for (j in seq_along(array2)){
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
      WithinGroupFactorsIC[1, 1:3, s] <- est_num(Xm,kmax = Rmax,type = array2[j])
      IC2[s] <- WithinGroupFactorsIC[1, 2, s]
  
      Llambda <- matrix(0, nrow = Nm, ncol = Smax)
      Lfhat <- matrix(0, nrow = T, ncol = Smax)
      LChat <- matrix(0, nrow = T, ncol = Nm)
  
      Sm <- min(15, ncol(vector))  # Προσαρμογή ώστε να μην υπερβαίνει το μήκος των δεδομένων
  
      if (Sm > 0) {
          Llambda[1:Nm, 1:Sm] <- sqrt(Nm) * Re(vector[, 1:Sm])
          Lfhat[1:T, 1:Sm] <- (1/Nm) * Xm %*% Llambda[1:Nm, 1:Sm]
          LChat[1:T, 1:Nm] <- Lfhat[1:T, 1:Sm] %*% t(Llambda[1:Nm, 1:Sm])
  }   else {
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
  
  
    }
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
out <- CP(list(X1, X2), rmax = 15, r0 = NULL,r = NULL, localfactor = TRUE, type = "ER")
Lambdahat <- as.matrix(out$loading_G[[1]])
globalCommon <- Fhat %*% t(Lambdahat)
Lambdahat <- as.matrix(out$loading_G[[2]])
globalCommon <- Fhat %*% t(Lambdahat)
Lambdahat <- as.matrix(do.call(rbind, out$loading_G))
globalCommon <- Fhat %*% t(Lambdahat)

plotHeatmapFull <- function(corr, main_title = "", colorbar = TRUE) {
  correlation <- corr
  dimMat <- ncol(correlation)
  color_palette <- colorRampPalette(c("white", "blue", "cyan", "yellow", "red"))(64)
  
  if (colorbar) {
    layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1.5))
  } else {
    layout(matrix(1))
  }
  
  par(mar = c(5, 5, 4, 2))
  
  image(
    x = 1:dimMat,
    y = 1:dimMat,
    z = t(correlation[dimMat:1, ]),
    col = color_palette,
    zlim = c(-1, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main_title
  )
  
  ticks <- seq(0, 88, by = 22)
  axis(1, at = ticks, labels = ticks)
  axis(2, at = ticks, labels = rev(ticks), las = 2)
  
  if (colorbar) {
    par(mar = c(5, 2, 4, 4))
    color_levels <- seq(-1, 1, length.out = 64)
    image(
      1,
      color_levels,
      matrix(color_levels, nrow = 1),
      col = color_palette,
      axes = FALSE,
      xlab = "",
      ylab = ""
    )
    axis(4, at = seq(-1, 1, by = 0.5), labels = seq(-1, 1, by = 0.5), las = 2)
  }
}


# Δημιουργία heatmap για τα κανονικοποιημένα δεδομένα
corr_data <- cor(X)
dev.new()
plotHeatmapFull(corr_data, main_title = "(a) Data Correlation", colorbar = TRUE)
Y_global <- X - globalCommon
corr_global <- cor(Y_global)
# Heatmap μετά την αφαίρεση του παγκόσμιου παράγοντα
dev.new()
plotHeatmapFull(corr_global, main_title = "(b) Residuals after Global Factor Removal", colorbar = TRUE)

# Heatmap μετά την αφαίρεση παγκόσμιου και τοπικού παράγοντα (GrFA)
residual_local <- do.call(cbind, out$residual)
corr_globloc <- cor(residual_local)
dev.new()
plotHeatmapFull(corr_globloc, main_title = "(c) Residuals after Global & Local Factor Removal", colorbar = TRUE)



