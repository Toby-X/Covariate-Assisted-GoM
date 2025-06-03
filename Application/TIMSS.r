install.packages(setdiff(
    c(
        "EdSurvey", "data.table", "stringr", "tidyverse",
        "cdmTools", "ggtern"
    ),
    rownames(installed.packages())
))

source("../utils.r")
source("../cogom.r")
library(EdSurvey)
library(data.table)
library(stringr)
library(tidyverse)
library(cdmTools)
library(ggtern)

## Cleaned TIMSS 2011 data for USA, grade 8 of 56 selected questions
load("TIMSS11.RData")
R_math <- R[, 1:26]
X_math <- X[, c(1:3, 5, 7, 9, 11, 13)]
X_math <- apply(X.math, 2, scale)

# Determine the number of latent profiles
paK(R_math)

K <- 3

# Without Covariates
res_null <- gom.svd(R_math, K)

# Parameter Tuning for Covariate-Assisted GoM
alpha_seq <- seq(from = 0, to = 1, length = 100)
err.math <- sapply(alpha_seq, .recovery.cv, R_math, X_math, K)
ggplot() +
    geom_smooth(aes(alpha_seq, err.math)) +
    geom_line(aes(alpha_seq, err.math), linewidth = .75) +
    labs(x = "alpha", y = "Recovery Error") +
    theme_minimal() +
    theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 16)
    )

alpha <- 0.45

# Covariate-Assisted GoM
res_cov <- gom.cov(R_math, X_math, K, alpha)

# Compare the results
## Heatmap of Theta
Theta.df_null <- data.frame(res_null$Theta)
rownames(Theta.df_null) <- paste0(ques.info[1:26], 1:26)
colnames(Theta.df_null) <- paste0("Profile ", 1:K)
Theta.df.null2 <- Theta.df_null %>%
    # as_tibble() %>%
    rownames_to_column("Questions") %>%
    pivot_longer(-Questions, names_to = "Profile", values_to = "value") %>%
    mutate(
        Profile = factor(Profile),
        Questions = factor(Questions, levels = unique(Questions))
    )

Theta.df_cov <- data.frame(res_cov$Theta)
rownames(Theta.df_cov) <- paste0(ques.info[1:26], 1:26)
colnames(Theta.df_cov) <- paste0("Profile ", 1:K)
library(tidyverse)
Theta.df.cov2 <- Theta.df_cov %>%
    # as_tibble() %>%
    rownames_to_column("Questions") %>%
    pivot_longer(-Questions, names_to = "Profile", values_to = "value") %>%
    mutate(
        Profile = factor(Profile),
        Questions = factor(Questions, levels = unique(Questions))
    )

Theta.df_full <- rbind(Theta.df.null2, Theta.df.cov2)
Theta.df_full$Covariates <- as.factor(
    rep(c("Without Covariates", "With Covariates"), each = nrow(Theta.df2))
)
pTheta <- ggplot(Theta.df_full, aes(Profile, Questions)) +
    geom_raster(aes(fill = value)) +
    geom_text(aes(label = round(value, 3))) +
    xlab("") +
    ylab(NULL) +
    scale_fill_distiller(palette = "Reds", direction = 1) +
    theme(
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 10),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5),
        strip.text = element_text(size = 16),
        strip.background = element_blank()
    ) +
    facet_wrap(~Covariates)
pTheta + scale_y_discrete(labels = function(x) str_wrap(x, width = 10))

## Ternary Plot of Pi
pi.dat <- data.frame(res_cov$Pi)
pi.dat$HomeResource <- df.usa$bsbgher
pi.dat$BulliedSchool <- df.usa$bsbgsbs
pi.dat$Like <- df.usa$bsbgsls
pi.dat$Value <- df.usa$bsbgsvs
pi.dat$Confident <- df.usa$bsbgscs
pi.dat$Engage <- df.usa$bsbgesl
pi.dat$GoodAt <- df.usa$bsbs19h
pi.dat$Income <- df.usa$bcbg05c

pi.dat2 <- data.frame(res_null$Pi)
pi.dat2$HomeResource <- df.usa$bsbgher
pi.dat2$BulliedSchool <- df.usa$bsbgsbs
pi.dat2$Like <- df.usa$bsbgsls
pi.dat2$Value <- df.usa$bsbgsvs
pi.dat2$Confident <- df.usa$bsbgscs
pi.dat2$Engage <- df.usa$bsbgesl
pi.dat2$GoodAt <- df.usa$bsbs19h
pi.dat2$Income <- df.usa$bcbg05c

### Covariance Results
ggtern(data = pi.dat, aes(X1, X2, X3)) +
    geom_point(aes(col = HomeResource), size = 2) +
    theme_hideticks() +
    theme_hidelabels() +
    xlab("Profile 1") +
    ylab("Profile 2") +
    zlab("Profile 3") +
    scale_color_distiller(palette = "YlOrRd") +
    theme(
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = .5)
    ) +
    labs(color = "value")
ggtern(data = pi.dat, aes(X1, X2, X3)) +
    theme_hideticks() +
    theme_hidelabels() +
    xlab("Profile 1") +
    ylab("Profile 2") +
    zlab("Profile 3") +
    scale_color_distiller(palette = "YlOrRd") +
    labs(color = "value")
ggtern(data = pi.dat, aes(X1, X2, X3)) +
    geom_point(aes(col = Engage), size = 2) +
    theme_hideticks() +
    theme_hidelabels() +
    xlab("Profile 1") +
    ylab("Profile 2") +
    zlab("Profile 3") +
    scale_color_distiller(palette = "YlOrRd")


### Without Covariance Results
ggtern(data = pi.dat2, aes(X1, X2, X3)) +
    geom_point(aes(col = factor(GoodAt)), size = 2) +
    theme_hideticks() +
    theme_hidelabels() +
    xlab("Profile 1") +
    ylab("Profile 2") +
    zlab("Profile 3") +
    scale_color_manual(values = c("#F1766D", "#FCBB44", "#839DD1", "#7A70B5")) +
    labs(color = "value")
ggtern(data = pi.dat2, aes(X1, X2, X3)) +
    geom_point(aes(col = HomeResource), size = 2) +
    theme_hideticks() +
    theme_hidelabels() +
    xlab("Profile 1") +
    ylab("Profile 2") +
    zlab("Profile 3") +
    scale_color_distiller(palette = "YlOrRd") +
    theme(
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = .5)
    ) +
    labs(color = "value")
ggtern(data = pi.dat2, aes(X1, X2, X3)) +
    geom_point(aes(col = factor(GoodAt)), size = 2) +
    theme_hideticks() +
    theme_hidelabels() +
    xlab("Profile 1") +
    ylab("Profile 2") +
    zlab("Profile 3") +
    scale_color_manual(values = c(
        "#F1766D",
        "#FCBB44", "#839DD1", "#7A70B5"
    )) +
    labs(color = "value")
