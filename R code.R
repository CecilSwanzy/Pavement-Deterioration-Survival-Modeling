# STEP 1 — LOAD PACKAGES
library(readxl)
library(data.table)
library(survival)
library(survminer)

# STEP 2 — READ DATA
nacka <- read_excel("nacka.xlsx")
ostersund <- read_excel("ostersund.xlsx")

setDT(nacka)
setDT(ostersund)

# Add municipality label
nacka$municipality     <- "Nacka"
ostersund$municipality <- "Ostersund"

# Make k_id globally unique 
nacka$k_id     <- paste0("N_", nacka$k_id)
ostersund$k_id <- paste0("O_", ostersund$k_id)

# Keep only necessary columns and combine
data <- rbind(
  nacka[, c("k_id", "year_week", "iri_mean",
            "variance_mean", "rcount_mean", "municipality")],
  ostersund[, c("k_id", "year_week", "iri_mean",
                "variance_mean", "rcount_mean", "municipality")]
)

setDT(data)

cat("Raw rows:", nrow(data), "\n")


# STEP 3 — CLEAN TYPES
data[, iri_mean := as.numeric(iri_mean)]
data[, variance_mean := as.numeric(variance_mean)]
data[, rcount_mean := as.numeric(rcount_mean)]

# Remove invalid rows
data <- data[
  !is.na(k_id) &
    !is.na(year_week) &
    !is.na(iri_mean)
]


# STEP 4 — AGGREGATE TO SECTION–WEEK LEVEL
data <- data[
  ,
  .(
    iri_mean = mean(iri_mean, na.rm = TRUE),
    variance_mean = mean(variance_mean, na.rm = TRUE),
    rcount_mean = mean(rcount_mean, na.rm = TRUE)
  ),
  by = .(k_id, year_week, municipality)
]

# Force uniqueness
data <- unique(data, by = c("k_id", "year_week", "municipality"))

dup_check <- data[, .N, by = .(k_id, year_week)][, max(N)]
cat("Max duplicates per section-week:", dup_check, "\n")


# STEP 5 — CREATE CLEAN TIME INDEX
data <- data[grepl("-", year_week)]

data[, c("year", "week") := tstrsplit(year_week, "-", fixed = TRUE)]

data[, year := as.numeric(year)]
data[, week := as.numeric(week)]

data <- data[!is.na(year) & !is.na(week)]

# Continuous weekly time
data[, time_index := year * 52 + week]

setorder(data, k_id, time_index)


# STEP 6 — DEFINE DETERIORATION EVENT
tau <- 3

data[, lag_iri := shift(iri_mean), by = k_id]

data[, event := as.integer(
  !is.na(lag_iri) &
    iri_mean >= tau &
    lag_iri < tau
)]

cat("Total deterioration events:", sum(data$event), "\n")
cat("Missing events:", sum(is.na(data$event)), "\n")


# STEP 7 — CREATE COUNTING PROCESS STRUCTURE
data[, start := shift(time_index), by = k_id]
data[, stop  := time_index]

data_surv <- data[!is.na(start)]

if(any(data_surv$stop <= data_surv$start)) {
  stop("Invalid survival intervals detected.")
}

cat("Max follow-up time:",
    max(data_surv$stop) - min(data_surv$start), "\n")


# STEP 8 — KAPLAN–MEIER 
first_event_section <- data_surv[
  ,
  {
    ev_rows <- which(event == 1)
    
    if (length(ev_rows) > 0) {
      list(
        event_time = stop[min(ev_rows)] - min(start),
        event_indicator = 1L
      )
    } else {
      list(
        event_time = max(stop) - min(start),
        event_indicator = 0L
      )
    }
  },
  by = .(k_id, municipality)
]

km_fit <- survfit(
  Surv(event_time, event_indicator) ~ municipality,
  data = first_event_section
)

ggsurvplot(
  km_fit,
  conf.int = TRUE,
  risk.table = TRUE
)


# STEP 9 — PREPARE STATIC COVARIATES

# Merge static variables from each municipality

# Nacka static variables 
nacka_static <- nacka[, .(
  k_id,
  funktion,
  area,
  langd,
  slitlager,
  prio,
  buss,
  roadtype,
  functional,
  farthinder
)]

# Östersund static variables 
ostersund_static <- ostersund[, .(
  k_id,
  typ,
  area,
  langd,
  belaggning,
  functional,
  farthinder
)]

# Harmonize names
setnames(nacka_static,
         c("funktion","slitlager"),
         c("road_type","surface_type"))

setnames(ostersund_static,
         c("typ","belaggning"),
         c("road_type","surface_type"))

static_data <- rbindlist(list(
  nacka_static[, .(k_id, road_type, surface_type,
                   functional, area, langd, farthinder)],
  ostersund_static[, .(k_id, road_type, surface_type,
                       functional, area, langd, farthinder)]
), fill = TRUE)

static_data <- unique(static_data, by = "k_id")

# Merge with first-event dataset
weibull_data <- merge(first_event_section,
                      static_data,
                      by = "k_id",
                      all.x = TRUE)

# Add CV summaries
cv_summary <- data[, .(
  variance_mean = mean(variance_mean, na.rm = TRUE),
  r_count_mean  = mean(rcount_mean, na.rm = TRUE)
), by = k_id]

weibull_data <- merge(weibull_data,
                      cv_summary,
                      by = "k_id",
                      all.x = TRUE)

# STEP 9.5 — FIX NUMERIC VARIABLES 

cat("\n========== FIXING NUMERIC VARIABLES BEFORE MODELING ==========\n")

# Convert to numeric
cat("\nConverting variables to numeric...\n")
weibull_data[, langd := as.numeric(langd)]
weibull_data[, area := as.numeric(area)]

cat("langd - NAs after conversion:", sum(is.na(weibull_data$langd)), "\n")
cat("area - NAs after conversion:", sum(is.na(weibull_data$area)), "\n")

# Remove rows with NA in critical variables
weibull_data <- weibull_data[!is.na(langd) & !is.na(area)]

cat("\nData rows after removing NAs:", nrow(weibull_data), "\n")

# Scale variables to improve convergence
cat("\nScaling continuous variables...\n")
weibull_data[, langd := scale(langd)]
weibull_data[, area := scale(area)]
weibull_data[, variance_mean := scale(variance_mean)]
weibull_data[, r_count_mean := scale(r_count_mean)]

cat("\nVariable statistics after scaling:\n")
print(summary(weibull_data[, .(langd, area, variance_mean, r_count_mean)]))

# Convert categorical variables to factors
weibull_data[, municipality := as.factor(municipality)]
weibull_data[, functional := as.factor(functional)]
weibull_data[, road_type := as.factor(road_type)]
weibull_data[, surface_type := as.factor(surface_type)]

# STEP 9.6 — DIAGNOSTIC CHECKS

cat("\n========== DIAGNOSTIC CHECKS ==========\n")

cat("\nCategorical variable distributions:\n")
print(table(weibull_data$municipality))
print(table(weibull_data$functional, weibull_data$municipality))
print(table(weibull_data$road_type, weibull_data$municipality))
print(table(weibull_data$surface_type, weibull_data$municipality))

# Check events by category
event_summary <- weibull_data[, .(
  sections = .N,
  events = sum(event_indicator)
), by = .(municipality, functional)]
print(event_summary)

# Check for missing values
missing_check <- weibull_data[, .(
  missing_functional = sum(is.na(functional)),
  missing_road_type = sum(is.na(road_type)),
  missing_surface_type = sum(is.na(surface_type)),
  missing_langd = sum(is.na(langd)),
  missing_area = sum(is.na(area))
)]
print(missing_check)

# Check correlation among numeric predictors
cat("\nCorrelation matrix:\n")
print(cor(weibull_data[, .(langd, area, variance_mean, r_count_mean)], 
          use = "complete.obs"))


# STEP 9.65 — FIX QUASI-SEPARATION

cat("\n========== FIXING RARE CATEGORIES ==========\n")

min_sections <- 50  # threshold


collapse_rare_levels <- function(varname) {
  
  # Tabulate sections and events by level
  tab <- weibull_data[, .(
    sections = .N,
    events   = sum(event_indicator)
  ), by = get(varname)]
  
  setnames(tab, "get", "level")
  
  # Identify problematic levels
  rare_levels <- tab[
    sections < min_sections | events == 0,
    level
  ]
  
  cat("\nCollapsing levels in", varname, ":\n")
  print(rare_levels)
  
  # Create fixed variable
  newvar <- paste0(varname, "_fixed")
  
  weibull_data[, (newvar) := as.character(get(varname))]
  
  weibull_data[get(varname) %in% rare_levels,
               (newvar) := "Other"]
  
  weibull_data[, (newvar) := factor(get(newvar))]
}

# Apply to categorical predictors
collapse_rare_levels("road_type")
collapse_rare_levels("surface_type")
collapse_rare_levels("functional")

cat("\nNew distributions after collapsing:\n")
print(table(weibull_data$road_type_fixed))
print(table(weibull_data$surface_type_fixed))
print(table(weibull_data$functional_fixed))


# STEP 9.7 — VERIFY MODELS RUN 

cat("\n========== VERIFYING MODELS CONVERGE ==========\n")

# Test models incrementally
cat("\nTesting model 1: municipality only\n")
cox_test1 <- coxph(Surv(event_time, event_indicator) ~ municipality, 
                   data = weibull_data)
print("✓ Model 1 converged")

cat("\nTesting model 2: + functional\n")
cox_test2 <- coxph(Surv(event_time, event_indicator) ~ municipality + functional, 
                   data = weibull_data)
print("✓ Model 2 converged")

cat("\nTesting model 3: + road_type\n")
cox_test3 <- coxph(Surv(event_time, event_indicator) ~ municipality + functional + road_type, 
                   data = weibull_data)
print("✓ Model 3 converged")

cat("\nTesting model 4: + surface_type\n")
cox_test4 <- coxph(Surv(event_time, event_indicator) ~ municipality + functional + road_type + surface_type, 
                   data = weibull_data)
print("✓ Model 4 converged")

cat("\nTesting model 5: + langd\n")
cox_test5 <- coxph(Surv(event_time, event_indicator) ~ municipality + functional + road_type + surface_type + langd, 
                   data = weibull_data)
print("✓ Model 5 converged")

cat("\nTesting model 6: + area\n")
cox_test6 <- coxph(Surv(event_time, event_indicator) ~ municipality + functional + road_type + surface_type + langd + area, 
                   data = weibull_data)
print("✓ Model 6 converged")

cat("\n✓✓✓ ALL MODELS CONVERGED SUCCESSFULLY ✓✓✓\n")

# STEP 10 — BASELINE COX (NO CV)

cat("\n========== BASELINE COX MODEL ==========\n")

cox_baseline <- coxph(
  Surv(event_time, event_indicator) ~
    municipality +
    functional_fixed +
    road_type_fixed +
    surface_type_fixed +
    langd +
    area,
  data = weibull_data
)

summary(cox_baseline)


# STEP 11 — EXTENDED COX (WITH CV)

cat("\n========== EXTENDED COX MODEL ==========\n")

cox_extended <- coxph(
  Surv(event_time, event_indicator) ~
    municipality +
    functional_fixed +
    road_type_fixed +
    surface_type_fixed +
    langd +
    area +
    variance_mean +
    r_count_mean,
  data = weibull_data
)

summary(cox_extended)

cat("\nModel comparison:\n")
anova(cox_baseline, cox_extended)


# COX PROPORTIONAL HAZARDS TEST

# Baseline model
ph_baseline <- cox.zph(cox_baseline)
print(ph_baseline)

# Extended model
ph_extended <- cox.zph(cox_extended)
print(ph_extended)


# PROPORTIONAL HAZARDS — BASELINE MODEL

# Remove GLOBAL row
vars_baseline <- rownames(ph_baseline$table)[-nrow(ph_baseline$table)]

n_vars <- length(vars_baseline)

# Layout: automatically choose rows/columns
n_col <- 2
n_row <- ceiling(n_vars / n_col)

par(mfrow = c(n_row, n_col),
    mar = c(5, 5, 3, 2),
    cex.lab = 1.3,
    cex.axis = 1.2,
    cex.main = 1.4)

for (v in vars_baseline) {
  
  plot(ph_baseline[v],
       main = v,
       xlab = "Time",
       ylab = "Scaled Schoenfeld Residuals",
       lwd = 2)
  
  abline(h = 0, col = "red", lwd = 2, lty = 2)
}

par(mfrow = c(1,1))


# PROPORTIONAL HAZARDS — EXTENDED MODEL

vars_extended <- rownames(ph_extended$table)[-nrow(ph_extended$table)]

n_vars_ext <- length(vars_extended)

n_col <- 2
n_row <- ceiling(n_vars_ext / n_col)

par(mfrow = c(n_row, n_col),
    mar = c(5, 5, 3, 2),
    cex.lab = 1.3,
    cex.axis = 1.2,
    cex.main = 1.4)

for (v in vars_extended) {
  
  plot(ph_extended[v],
       main = v,
       xlab = "Time",
       ylab = "Scaled Schoenfeld Residuals",
       lwd = 2)
  
  abline(h = 0, col = "red", lwd = 2, lty = 2)
}

par(mfrow = c(1,1))


# HANDLING THE PH ASSUMPTION VIOLATION

# Add time-varying effects

cox_extended_tv <- coxph(
  Surv(event_time, event_indicator) ~
    functional_fixed +
    road_type_fixed +
    surface_type_fixed +
    langd +
    area +
    variance_mean +
    r_count_mean +
    tt(variance_mean) +
    tt(r_count_mean) +
    strata(municipality),
  data = weibull_data,
  tt = function(x, t, ...) x * log(t)
)



# Stratify by municipality

cox_extended_fixed <- coxph(
  Surv(event_time, event_indicator) ~
    functional_fixed +
    road_type_fixed +
    surface_type_fixed +
    langd +
    area +
    variance_mean +
    r_count_mean +
    strata(municipality),
  data = weibull_data
)

# STEP 2 — Re-test PH

ph_fixed <- cox.zph(cox_extended_fixed)
print(ph_fixed)


# STEP 12 — WEIBULL MODELS

cat("\n========== WEIBULL MODELS ==========\n")

weibull_baseline <- survreg(
  Surv(event_time, event_indicator) ~
    municipality +
    functional_fixed +
    road_type_fixed +
    surface_type_fixed +
    langd +
    area,
  data = weibull_data,
  dist = "weibull"
)

weibull_extended <- survreg(
  Surv(event_time, event_indicator) ~
    municipality +
    functional_fixed +
    road_type_fixed +
    surface_type_fixed +
    langd +
    area +
    variance_mean +
    r_count_mean,
  data = weibull_data,
  dist = "weibull"
)

cat("\nAIC comparison:\n")
print(AIC(weibull_baseline, weibull_extended))

summary(weibull_extended)

# STEP 13 — BOOTSTRAP UNCERTAINTY 

set.seed(123)

B <- 300
boot_medians <- numeric(B)

section_ids <- unique(weibull_data$k_id)

for (b in 1:B) {
  
  # Resample sections with replacement
  sampled_ids <- sample(section_ids,
                        size = length(section_ids),
                        replace = TRUE)
  
  boot_data <- weibull_data[k_id %in% sampled_ids]
  
  # Refit Weibull model
  boot_model <- tryCatch(
    survreg(
      Surv(event_time, event_indicator) ~ municipality +
        functional_fixed +
        road_type_fixed +
        surface_type_fixed +
        langd +
        area +
        variance_mean +
        r_count_mean,
      data = boot_data,
      dist = "weibull"
    ),
    error = function(e) return(NULL)
  )
  
  if (!is.null(boot_model)) {
    
    shape_b <- 1 / boot_model$scale
    
    #  SAFE PREDICTION APPROACH 
    # Predict linear predictor for each section
    eta_hat <- predict(boot_model, type = "lp")
    
    # Compute median survival per section
    median_pred <- exp(eta_hat) * (log(2))^(1/shape_b)
    
    # Store overall median predicted maintenance time
    boot_medians[b] <- median(median_pred)
    
  } else {
    boot_medians[b] <- NA
  }
}

# Remove failed iterations
valid_medians <- boot_medians[!is.na(boot_medians)]

cat("Number of successful iterations:", length(valid_medians), "\n")
cat("Success rate:",
    round(length(valid_medians)/B * 100, 1), "%\n")

# Confidence interval
ci_boot <- quantile(valid_medians, c(0.025, 0.975))

cat("Bootstrap 95% CI for median maintenance time:\n")
print(ci_boot)

# Distribution plot
hist(valid_medians,
     breaks = 30,
     main = "Bootstrap Distribution of Median Maintenance Time",
     xlab = "Predicted Median Time (weeks)")

abline(v = median(valid_medians), lwd = 2)

# Normality diagnostic 
qqnorm(valid_medians)
qqline(valid_medians)

length(section_ids)
summary(valid_medians)
sd(valid_medians)


# FAST COX BOOTSTRAP UNCERTAINTY 

set.seed(123)

B <- 200
boot_medians_cox <- numeric(B)

section_ids <- unique(weibull_data$k_id)

for (b in 1:B) {
  
  sampled_ids <- sample(section_ids,
                        size = length(section_ids),
                        replace = TRUE)
  
  boot_data <- weibull_data[k_id %in% sampled_ids]
  
  boot_data <- boot_data[complete.cases(boot_data[, c(
    "event_time",
    "event_indicator",
    "municipality",
    "functional_fixed",
    "road_type_fixed",
    "surface_type_fixed",
    "langd",
    "area",
    "variance_mean",
    "r_count_mean"
  )]), ]
  
  if (nrow(boot_data) == 0) {
    boot_medians_cox[b] <- NA
    next
  }
  
  boot_model <- tryCatch(
    coxph(
      Surv(event_time, event_indicator) ~ municipality +
        functional_fixed +
        road_type_fixed +
        surface_type_fixed +
        langd +
        area +
        variance_mean +
        r_count_mean,
      data = boot_data
    ),
    error = function(e) return(NULL)
  )
  
  if (!is.null(boot_model)) {
    
    # Vectorized survival computation
    sf_all <- tryCatch(
      survfit(boot_model, newdata = boot_data),
      error = function(e) return(NULL)
    )
    
    if (!is.null(sf_all)) {
      
      medians_section <- numeric(nrow(boot_data))
      
      # sf_all$surv is matrix: rows=time, cols=individuals
      for (i in 1:ncol(sf_all$surv)) {
        
        surv_i <- sf_all$surv[, i]
        time_i <- sf_all$time
        
        idx <- which(surv_i <= 0.5)[1]
        
        if (!is.na(idx)) {
          medians_section[i] <- time_i[idx]
        } else {
          medians_section[i] <- NA
        }
      }
      
      boot_medians_cox[b] <- median(medians_section, na.rm = TRUE)
      
    } else {
      boot_medians_cox[b] <- NA
    }
    
  } else {
    boot_medians_cox[b] <- NA
  }
}

valid_medians_cox <- boot_medians_cox[!is.na(boot_medians_cox)]

cat("Number of successful iterations:", length(valid_medians_cox), "\n")
cat("Success rate:",
    round(length(valid_medians_cox)/B * 100, 1), "%\n")

ci_boot_cox <- quantile(valid_medians_cox, c(0.025, 0.975))

cat("Bootstrap 95% CI (Cox):\n")
print(ci_boot_cox)

summary(valid_medians_cox)
sd(valid_medians_cox)

# Cox Survival Plot
sf_check <- survfit(coxph(
  Surv(event_time, event_indicator) ~ municipality +
    functional_fixed +
    road_type_fixed +
    surface_type_fixed +
    langd +
    area +
    variance_mean +
    r_count_mean,
  data = weibull_data
))

# Create plot with proper axis labels
plot(sf_check, 
     xlab = "Time", 
     ylab = "Survival Probability",
     main = "Cox Proportional Hazards Survival Curve")

# Add horizontal line at 0.5
abline(h = 0.5, col = "red")

# Add legend
legend("topright", 
       legend = c("Survival Curve", "Median Survival (0.5)"),
       col = c("black", "red"),
       lty = c(1, 1),
       lwd = c(1, 1))

# STEP 14 — SUMMARY OUTPUT

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Final dataset dimensions:", dim(weibull_data), "\n")
cat("Number of sections:", length(unique(weibull_data$k_id)), "\n")
cat("Total events:", sum(weibull_data$event_indicator), "\n")
cat("Censored:", sum(weibull_data$event_indicator == 0), "\n")


levels(weibull_data$road_type_fixed)
levels(weibull_data$functional_fixed)

contrasts(weibull_data$road_type_fixed)
summary(cox_baseline)


