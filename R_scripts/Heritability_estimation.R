##############
#Using SOMMER#
##############
library("bigmemory")
library("sommer")
library("biganalytics")
library(tidyverse)


data = read.table("prt_lines.txt",header = T, sep = " ")
data = read.csv("final_MD.csv", header = TRUE)

grm <- as.matrix(read.table("md_data.rel", header = FALSE))
ids <- read.table("md_data.rel.id")
clean_ids <- sub("_(.*)", "", ids$V2)
rownames(grm) <- clean_ids
colnames(grm) <- clean_ids

#Remove missing values
data = data %>% filter(!is.na(MD))

## Single environment 
## To calculate H and h
ans1 <- mmer(BBD~Loc_Year, 
             random= ~ LINE + LINE:Loc_Year, data=data, verbose = FALSE) #Univariate homogeneous variance model


summary(ans1)$varcomp

# Variance components from the model output

#replace with your values
var_line <- 83.44914        # Genetic variance (Line)
var_gxe <- 0.00000      # Genotype by environment interaction variance (Line:Env)
var_residual <- 33.63771

##Broad sense heritability##

# Total phenotypic variance
var_phenotypic <- var_line + var_gxe + var_residual

# Broad-sense heritability (H^2)
heritability <- (var_line + var_gxe) / var_phenotypic

# Output the heritability
heritability


## Narrow sense heritability without GRM##

# Total phenotypic variance
var_phenotypic <- var_line + var_gxe + var_residual

# Narrow-sense heritability (h^2)
h2 <- var_line / var_phenotypic

# Output the narrow-sense heritability
h2

## Narrow sense heritability with GRM##
#Multi-environment (doesn't work properly for mine)

ans2r <- mmer(BBD~Loc_Year, random= ~LINE + vsr(dsr(Loc_Year),LINE), rcov= ~ vsr(dsr(Loc_Year),units), data=data, verbose = FALSE)

ans2g <- mmer(
  PRT ~ Loc_Year, 
  random = ~ vs(LINE, Gu = A) + vsr(dsr(Loc_Year), LINE),  # Line uses GRM, GxE is still modeled
  rcov   = ~ vsr(dsr(Loc_Year), units),
  data   = data,
  verbose = FALSE
)

## To calculate hg
ans2g <- mmer(
  MD ~ Loc_Year,
  random = ~ vsr(LINE, Gu=grm),  # grm = additive relationship matrix
  rcov = ~ units,
  data = data
)

summary(ans2g)$varcomp

# Variance components
var_additive <- 71.58283 #LINE
var_residual <- 58.11967 #units

# Total phenotypic variance
var_total <- var_additive + var_residual

# Narrow-sense heritability
h2 <- var_additive / var_total

# Print result
h2

##################
####Using lme4####
##################
library(lme4)

data = read.table("prt_lines.txt",header = T, sep = ",")

#Remove missing values
data = data %>% filter(!is.na(PRT))

#make a new variable location_year
data$ENV = paste0(data$Location,"_",data$Year)

#Run model

model <- lmer(PRT ~ ENV + ENV:Line + (1 | Line), data = data) # with GxE interaction as fixed effect
model <- lmer(PRT ~ ENV + (1 | Line), data = data)
model <- lmer(PRT ~ ENV + (1 | Line) + (1 | Line:ENV), data = data) # with GxE interaction as random effect same a sommer effect BUT IT DOESN'T WORK  as the data is sparse

# Print the model summary to view the variance components
summary(model)

# Extract the variance components
varcomp <- as.data.frame(VarCorr(model))

## Broad sense heritability ##

# Calculate genetic variance (VG) which includes variance due to Line and Location
VG <- varcomp[1, "vcov"] + varcomp[2, "vcov"]  # Genotypic variance (Line + Location)

# Extract the residual (error) variance (Ve)
Ve <- attr(VarCorr(model), "sc")^2  # Residual variance

# Calculate the total phenotypic variance (VP)
VP <- VG + Ve

# Calculate broad-sense heritability (H2)
H2 <- VG / VP

# Output the broad-sense heritability
H2


## narrow sense heritability ##


# Additive genetic variance (from Line)
additive_variance <- varcomp[varcomp$grp == "Line", "vcov"]

# Residual variance (from the model's scale)
residual_variance <- attr(VarCorr(model), "sc")^2  # Residual variance

# Phenotypic variance is the sum of genetic variance and residual variance
phenotypic_variance <- additive_variance + residual_variance

# Calculate narrow-sense heritability (h^2)
h2 <- additive_variance / phenotypic_variance

# Print narrow-sense heritability
print(paste("Narrow-sense heritability: ", round(h2, 3)))

