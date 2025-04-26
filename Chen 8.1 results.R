library(GrFA)
library(gridExtra)
# Διαγράφει όλες τις μεταβλητές από το περιβάλλον εργασίας.
# Κλείνει όλα τα ανοιχτά γραφήματα.
graphics.off(); rm(list = ls(all = TRUE)); cat("\014"); getwd()
# Εισαγωγή δεδομένων
csv_in_developed <- "dev1.csv"
data_developed <- read.csv(csv_in_developed,header=F)  # Διαβάζει τα δεδομένα από το αρχείο CSV για την ανεπτυγμένη ομάδα.
Ncol1 <- ncol(data_developed)  # Αποθηκεύει τον αριθμό των στηλών για την ανεπτυγμένη ομάδα.

csv_in_emerging <- "eme1.csv"
data_emerging <- read.csv(csv_in_emerging, header=F)  # Διαβάζει τα δεδομένα από το αρχείο CSV για την αναπτυσσόμενη ομάδα.
Ncol2 <- ncol(data_emerging)  # Αποθηκεύει τον αριθμό των στηλών για την αναπτυσσόμενη ομάδα.

# Κανονικοποίηση δεδομένων
X0 <- cbind(data_developed, data_emerging)  # Συνδυάζει τα δεδομένα από τις δύο ομάδες.
X <- scale(X0, center = TRUE, scale = TRUE)  # Κανονικοποιεί τις στήλες (μέσος όρος 0, τυπική απόκλιση 1).
mean(apply(X, 2, sd))  # Ελέγχει αν η τυπική απόκλιση όλων των στηλών είναι 1.
sum(colMeans(X))  # Ελέγχει αν ο μέσος όρος όλων των στηλών είναι κοντά στο 0.
X1=X[,1:37]
X2=X[,38:88]





# Υπολογισμός PCA για τοπικούς παράγοντες
pca_results <- prcomp(X, center = TRUE, scale. = TRUE)

# Εμφάνιση ιδιοτιμών (ποσοστά διακύμανσης)
pca_var <- pca_results$sdev^2
explained_variance <- pca_var / sum(pca_var)

cat("Ποσοστά εξηγούμενης διακύμανσης:\n")
print(explained_variance)

# Επιλογή αριθμού παραγόντων με βάση την εξήγηση του 90% της διακύμανσης
cumulative_variance <- cumsum(explained_variance)
num_factors <- which(cumulative_variance >= 0.90)[1]

cat("Αριθμός παραγόντων που εξηγούν το 90% της διακύμανσης:", num_factors, "\n")

# Τοπικοί παράγοντες
local_factors <- pca_results$x[, 1:num_factors]

# Εμφάνιση πρώτων τοπικών παραγόντων
cat("Πρώτοι τοπικοί παράγοντες:\n")
head(local_factors)

# Υπολογισμός παγκόσμιων παραγόντων με Circular Projection
# Υπολογισμός πίνακα συνδιακύμανσης τοπικών παραγόντων
cov_matrix_local <- cov(local_factors)

# Υπολογισμός ιδιοτιμών και ιδιοδιανυσμάτων για παγκόσμιους παράγοντες
eig_results_global <- eigen(cov_matrix_local)

global_eigenvalues <- eig_results_global$values
global_eigenvectors <- eig_results_global$vectors

cat("Ιδιοτιμές παγκόσμιων παραγόντων:\n")
print(global_eigenvalues)

# Επιλογή παγκόσμιων παραγόντων (με βάση τις 5 πρώτες ιδιοτιμές όπως αναφέρεται)
global_factors <- local_factors %*% global_eigenvectors[, 1:5]

# Εμφάνιση πρώτων παγκόσμιων παραγόντων
cat("Πρώτοι παγκόσμιοι παράγοντες:\n")
head(global_factors)

X3 <- as.matrix(X1)
X4 <- as.matrix(X2)
anyNA(X3)
anyNA(X4)
i_values <- 1:floor(Ncol1/2)

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
i_values <- 1:floor(Ncol2/2)

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
i_values <- 1:floor(Ncol1/2)

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
i_values <- 1:floor(Ncol2/2)

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


library(ggplot2)
library(tidyr)
library(dplyr)
T <- nrow(X)                # Αριθμός γραμμών (χρονικές στιγμές).
N1 <- Ncol1                 # Αριθμός στηλών για την ανεπτυγμένη ομάδα.
N2 <- Ncol2                 # Αριθμός στηλών για την αναπτυσσόμενη ομάδα.
N <- c(N1, N2)              # Πίνακας που περιέχει τον αριθμό των στηλών για κάθε ομάδα.
M <- 2
# Ορισμός πινάκων προβολής
Proj <- array(0, dim = c(T, T, M))  # Δημιουργία τρισδιάστατου πίνακα για αποθήκευση προβολών κάθε ομάδας.
CirProjection <- diag(T)           # Αρχικός κυκλικός πίνακας προβολής (μοναδιαίος πίνακας).

# Προετοιμασία μεταβλητών για υπολογισμό παραγόντων
WithinGroupFactorsIC <- matrix(0, nrow = M, ncol = 1)  # Πίνακας για αποθήκευση αποτελεσμάτων για κάθε ομάδα.
kkk <- 0  # Μεταβλητή για την παρακολούθηση των στηλών.

# Υπολογισμός παραγόντων για κάθε ομάδα
for (s in 1:M) {
  Nm <- N[s]                               # Αριθμός στηλών για την ομάδα `s`.
  Xm <- X[, (kkk + 1):(kkk + Nm)]          # Υποσύνολο του X για την ομάδα `s`.
  
  eig_result <- eigen(crossprod(Xm) / (Nm * T))  # Υπολογισμός ιδιοτιμών και ιδιοδιανυσμάτων.
  eigenvalues <- eig_result$values              # Αποθήκευση ιδιοτιμών.
  
  position <- which.max(eigenvalues[1:(length(eigenvalues) - 1)] / 
                          eigenvalues[2:length(eigenvalues)])  # Εύρεση της θέσης της μέγιστης αναλογίας ιδιοτιμών.
  WithinGroupFactorsIC[s] <- position           # Αποθήκευση της θέσης για την ομάδα `s`.
  
  # Υπολογισμός συνιστωσών
  Sm <- min(15, Nm)                            # Ανώτατος αριθμός συνιστωσών (μέγιστο 15 ή Nm).
  Llambda <- sqrt(Nm) * eig_result$vectors[, 1:Sm]  # Υπολογισμός φορτίων παραγόντων.
  Lfhat <- (1 / Nm) * Xm %*% Llambda               # Υπολογισμός εκτιμήσεων παραγόντων.
  LChat <- Lfhat %*% t(Llambda)                   # Ανακατασκευή των δεδομένων.
  Ldevi <- norm(LChat - Xm, type = "F") / sqrt(Nm * T)  # Απόκλιση της ανακατασκευής.
  
  KmHat <- Lfhat  # Αποθήκευση των εκτιμήσεων παραγόντων.
  if (s == 1) {
    K1hat <- KmHat  # Αποθήκευση για την πρώτη ομάδα.
  }
  
  Proj[, , s] <- KmHat %*% solve(t(KmHat) %*% KmHat) %*% t(KmHat)  # Δημιουργία πίνακα προβολής για την ομάδα.
  CirProjection <- CirProjection %*% Proj[, , s]  # Ενημέρωση κυκλικού πίνακα προβολής.
  kkk <- kkk + Nm  # Ενημέρωση της εκκίνησης για τις στήλες.
}

# Υπολογισμός ιδιοδιανυσμάτων
eig_CirProjection <- eigen(CirProjection)
global <- as.data.frame(-Re(eig_CirProjection$vectors[, 1]))  # Μόνο το πραγματικό μέρος
colnames(global) <- paste0("V", 1)  # Ονοματοδοσία στηλών

# Προσθήκη άξονα χρόνου
global$Time <- 1:nrow(global)

# Μετατροπή σε long format για ggplot
global_long <- pivot_longer(global, cols = starts_with("V"), names_to = "Factor", values_to = "Value")

# Δημιουργία γραφήματος
library(ggplot2)
library(dplyr)
library(tidyr)

# Ορισμός των breaks και των labels για τα τρίμηνα
my_breaks <- c(1, 12, 28, 44, 60, 76, 92, 108)   # Αντιστοιχούν στις θέσεις (indices) των τριμήνων
my_labels <- c("Q3/98", "Q2/01", "Q1/04", "Q4/06", "Q3/09", "Q2/12", "Q1/15", "Q4/17")

plot_global <- ggplot(global_long, aes(x = Time, y = Value)) +
  geom_line(color = "black", size = 0.8) +
  facet_wrap(~ Factor, ncol = 2, scales = "free_y") +
  labs(
    title = "Σχήμα 1(a):  Πρώτος Παγκόσμιος Παράγοντας",
    x = "Χρόνος", 
    y = "Τιμή Παράγοντα"
  ) +
  scale_x_continuous(
    limits = c(1, 108),         # Ορίζει τα όρια ώστε να συμπεριληφθούν όλα τα τρίμηνα
    breaks = my_breaks, 
    labels = my_labels,
    expand = c(0.01, 0.01)       # Μειώνει τα κενά πριν/μετά τον πρώτο/τελευταίο δείκτη
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),   # Έντονοι τίτλοι των facet
    panel.grid.minor = element_blank(),         # Αφαιρεί το λεπτό πλέγμα
    axis.text.x = element_text(size = 7)          # Μικρότερο μέγεθος γραμματοσειράς
    # Εναλλακτικά: axis.text.x = element_text(size = 7, angle = 45, hjust = 1)
  )

print(plot_global)
ACP <- function(y, reference_group, rmax = 8, r0 = NULL, r = NULL, r_aug = NULL, localfactor = FALSE, type = "IC3") {
  M <- length(y)
  T <- nrow(y[[1]])
  Nm <- sapply(y, ncol)
  K <- list()
  Proj_Mat <- list()
  
  for (m in 1:M) {
    K[[m]] <- FA(y[[m]], r = rmax)$F
    Proj_Mat[[m]] <- K[[m]] %*% solve(t(K[[m]]) %*% K[[m]]) %*% t(K[[m]])
  }
  
  Aug_Proj_Mat <- Proj_Mat[[reference_group]]
  for (m in setdiff(1:M, reference_group)) {
    Aug_Proj_Mat <- Aug_Proj_Mat %*% Proj_Mat[[m]]
  }
  
  Aug_Proj_Mat <- t(Aug_Proj_Mat) %*% Aug_Proj_Mat
  eig_deco <- eigen(Aug_Proj_Mat)
  eig_vec <- eig_deco$vectors[, 1:rmax]
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
    if (!(r0 %% 1 == 0) | r0 <= 0) {
      stop("invalid 'r0' input")
    }
    r0hat <- r0
  }
  
  Ghat <- as.matrix(eig_deco$vectors[, 1:r0hat]) * sqrt(T)
  Proj_G <- Ghat %*% solve(t(Ghat) %*% Ghat) %*% t(Ghat)
  loading_G <- list()
  for (m in 1:M) {
    loading_G[[m]] <- 1/T * t(y[[m]]) %*% Ghat
  }
  
  # Εκτίμηση των Augmented Global Factors
  eig_aug <- eigen(Aug_Proj_Mat)
  
  if (is.null(r_aug)) {
    r_aug <- sum(eig_aug$values > 1e-4)  # Αν δεν έχει δοθεί, υπολογίζεται από τα σημαντικά eigenvalues
  } else {
    if (!(r_aug %% 1 == 0) | r_aug <= 0) {
      stop("invalid 'r_aug' input")
    }
  }
  
  AGhat <- as.matrix(eig_aug$vectors[, 1:r_aug]) * sqrt(T)  # Augmented Global Factors matrix
  
  rhat <- rep(0, M)  # Αρχικοποίηση του αριθμού των τοπικών παραγόντων για κάθε group
  
  if (localfactor == FALSE) {
    res <- list(r0hat = r0hat, rhat = rhat, r_aug = r_aug, rho = rho, Ghat = Ghat, AGhat = AGhat, loading_G = loading_G)
  } else {
    Fhat <- list()
    loading_F <- list()
    y_proj_G <- lapply(y, function(x) x - Proj_G %*% x)
    
    if (is.null(r)) {
      for (m in 1:M) {
        rhat[m] <- est_num(y_proj_G[[m]], kmax = rmax - r0hat, type = type)
        fit <- FA(y_proj_G[[m]], r = rhat[m])
        Fhat[[m]] <- fit$F
        loading_F[[m]] <- fit$L
      }
    } else {
      if (!(all(r %% 1 == 0) && all(r >= 0))) {
        stop("invalid r input")
      }
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
      r0hat = r0hat,  # Αριθμός παγκόσμιων παραγόντων
      rhat = rhat,    # Αριθμός τοπικών παραγόντων ανά group
      r_aug = r_aug,  # Αριθμός Augmented Global Factors
      rho = rho, 
      Ghat = Ghat, 
      AGhat = AGhat,  # Augmented Global Factors matrix
      Fhat = Fhat, 
      loading_G = loading_G, 
      loading_F = loading_F, 
      residual = e
    )
  }
  
  class(res) <- "GFA"
  return(res)
}
out2 <- ACP(list(X1,X2),reference_group = 1,rmax = 15,r0 = 5,r = rep(1,2),r_aug = 5,localfactor = TRUE,type = "ER")
# Εκτύπωση στο δεύτερο panel
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(purrr)

# Καθορισμός των σωστών breakpoints και labels ώστε να περιλαμβάνει και το Q4/17
my_breaks <- seq(1, nrow(out2$AGhat), length.out = 8)  # 8 σημεία
my_labels <- c("Q3/98", "Q2/01", "Q1/04", "Q4/06", "Q3/09", "Q2/12", "Q1/15", "Q4/17")  # 8 ετικέτες

# ---- Προετοιμασία των 5 πρώτων Augmented Global Factors ----
df_aug <- data.frame(Time = seq_len(nrow(out2$AGhat)), as.data.frame(-out2$AGhat[, 1:5]))
colnames(df_aug) <- c("Time", paste0("Aug_Global_Factor_", 1:5))

df_aug_long <- pivot_longer(df_aug, cols = -Time, names_to = "Factor", values_to = "Value")

# ---- Προετοιμασία του Global Factor ----


# ---- Γράφημα των 5 πρώτων Augmented Global Factors ----
plot_aug <- ggplot(df_aug_long, aes(x = Time, y = Value, color = Factor)) +
  geom_line(size = 1) +
  labs(title = "Augmented Global Factors (Top 5)", x = "Χρόνος", y = "Τιμή Παράγοντα") +
  scale_x_continuous(breaks = my_breaks, labels = my_labels, expand = c(0.01, 0.01)) +
  theme_minimal() +
  theme(
    legend.title = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1)
  )

print(plot_aug)  # Εκτύπωση στο πρώτο panel

# ---- Γράφημα του Global Factor και του 1ου Augmented Global Factor ----


library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Καθορισμός των σωστών breakpoints και labels ώστε να περιλαμβάνει και το Q4/17
my_breaks <- seq(1, nrow(out2$AGhat), length.out = 8)  # 8 σημεία
my_labels <- c("Q3/98", "Q2/01", "Q1/04", "Q4/06", "Q3/09", "Q2/12", "Q1/15", "Q4/17")  # 8 ετικέτες

# ---- Προετοιμασία των 5 πρώτων Augmented Global Factors ----
df_aug <- data.frame(Time = seq_len(nrow(-out2$AGhat)), as.data.frame(-out2$AGhat[, 1:5]))
colnames(df_aug) <- c("Time", paste0("AG", 1:5))  # Μετονομασία σε V1, V2, ..., V5

df_aug_long <- pivot_longer(df_aug, cols = -Time, names_to = "Factor", values_to = "Value")

# ---- Δημιουργία λίστας με τα 5 γραφήματα ----
plot_list <- list()

for (i in 1:5) {
  factor_name <- paste0("AG", i)
  
  plot_list[[i]] <- ggplot(filter(df_aug_long, Factor == factor_name), aes(x = Time, y = Value)) +
    geom_line(color = "black", size = 0.8) +  # Μαύρο χρώμα για όλες τις γραμμές
    labs(title = factor_name, x = NULL, y = "Τιμή Παράγοντα") +
    scale_x_continuous(breaks = my_breaks, labels = my_labels, expand = c(0.01, 0.01)) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)  # Κεντραρισμένος τίτλος
    )
}

# ---- Εμφάνιση των 5 γραφημάτων σε διάταξη 3x2 ----
grid.arrange(
  plot_list[[1]], plot_list[[2]],
  plot_list[[3]], plot_list[[4]],
  plot_list[[5]], nrow = 3, ncol = 2, bottom = "Χρόνος"
)




library(GrFA)
library(gridExtra)

# --------------------------------------------
# Βήμα 1: Φόρτωση και προετοιμασία δεδομένων
# --------------------------------------------
# Emerging group

# --------------------------------------------
# Βήμα 2: Κλήση της συνάρτησης ACP
# --------------------------------------------
y <- list(X1, X2)
acp_result <- ACP(
  y = y,
  reference_group = 1,  # Developed ως reference
  rmax = 15,            # Μέγιστος αριθμός factors
  r_aug = 5             # Augmented factors
)
library(GrFA)
library(ggplot2)
# Αριθμός παραγόντων

# Αρχικός κυκλικός πίνακας



# Εξαγωγή αποτελεσμάτων
Ghat <- acp_result$Ghat[,1:5]     # Global Factors (Tn x R0)
AGhat <- Re(eig_CirProjection$vectors[,1:5])
# Augmented Global Factors (Tn x Raug)

# --------------------------------------------
# Βήμα 3: Βελτιστοποιημένες συναρτήσεις διακύμανσης
# --------------------------------------------
calculate_variance <- function(data, factors) {
  # Βεβαιωθείτε ότι οι παράγοντες έχουν σωστές διαστάσεις
  if (nrow(data) != nrow(factors)) stop("Μη συμβατές διαστάσεις")
  
  # Υπολογισμός φορτίων και προβολής
  loadings <- t(data) %*% factors %*% solve(t(factors) %*% factors)
  projected <- factors %*% t(loadings)
  
  # Υπολογισμός διακυμάνσεων
  total_var <- sum(apply(data, 2, var))
  explained_var <- sum(apply(projected, 2, var))
  explained_var / total_var
}

# --------------------------------------------
# Βήμα 4: Υπολογισμός τελικών τιμών
# --------------------------------------------
# Within-Group Variance
within_dev <- summary(prcomp(X1))$importance[2, 1] |> round(3) # 0.445
within_eme <- summary(prcomp(X2))$importance[2, 1] |> round(3) # 0.256

# Global & Augmented Variance (χωρίς προσαρμογή)
global_dev <- calculate_variance(X1, Ghat) |> round(3)
global_eme <- calculate_variance(X2, Ghat) |> round(3)
aug_dev <- calculate_variance(X1, AGhat) |> round(3)
aug_eme <- calculate_variance(X2, AGhat) |> round(3)

# --------------------------------------------
# Βήμα 5: Δημιουργία πίνακα
# --------------------------------------------
variance_table <- data.frame(
  Group = c("Developed (group 1):", "Emerging (group 2):"),
  First_Within_Group = c(within_dev, within_eme),
  Global = c(global_dev, global_eme),
  Aug_Global = c(aug_dev, aug_eme)
)
table_grob <- gridExtra::tableGrob(
  variance_table,
  rows = NULL,
  cols = c("", "The first within-group", "Global", "Aug. global"),
  theme = gridExtra::ttheme_minimal(
    core = list(bg_params = list(fill = c("grey95", "white"))),
    colhead = list(bg_params = list(fill = "grey80"))
  ))
# --------------------------------------------
# Βήμα 6: Εμφάνιση πίνακα
#grid.arrange(
library(ggplot2)
ggplot() +
  theme_void() +
  annotation_custom(table_grob) +
  labs(title = "Table 4. Variances explained by the first within-group factor, the five global factors and the five augmented global factors.")

write.csv(out2$AGhat[,1:5], "globalFactorsAug.csv", row.names = FALSE)
write.csv(Re(eig_CirProjection$vectors[,1:5]), "globalFactors.csv", row.names = FALSE)