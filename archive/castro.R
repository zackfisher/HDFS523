# SIMULATE DATASETS A AND B
# based on RM ANOVA specification model (van Breukelen, 2013)
# Set seed to allow replication
set.seed (1234)
# Create ID variable (N = 100)
id <- c(1:100, 1:100)
# Create a time variable (0 = time 1, 1 = time 2)
time <- c(rep(0, 100), rep(1, 100))
# Create grouping variable with equal number cohabiting (1) and not (0)
group <- rep(c(0,1), 100)
# Put all variables in a dataframe
mydf <- data.frame(id, group, time)

# Generate relationship satisfaction score
# To simulate dataset B, replace the value of -1 in the line below with 0
mydf$score = 4.5 + -1*mydf$group + 0*mydf$time + 
             0*mydf$group*mydf$time + rnorm(200, 0, 0.3)

# Make into wide format for computing difference score
mydf <- reshape(
  mydf, 
  v.names = 'score', 
  timevar = "time",
  idvar = "id", 
  direction= "wide"
)

# Assign appropriate variable names
names(mydf) = c('id', 'Cohabit', 'rs1', 'rs2')

# Compute difference variable
mydf$rsdiff = mydf$rs2 - mydf$rsl

# Declare cohabiting variable as a factor (categorical variable)
mydf$Cohabit = factor(mydf$Cohabit)

# VISUALIZE THE DATA
library(ggplot2)
ggplot(mydf, aes(x = rsl, y = rs2, shape= Cohabit)) + 
  geom_point(size = 3) +
  geom_smooth(method = lm, se = F) +
  xlab("Relationship Satisfaction Time 1") +
  ylab("Relationship Satisfaction Time 2")


# FIT MODELS
# ANCOVA model
mylml <- lm(rs2 ~ rsl + Cohabit, data = mydf)
summary (mylml)


# Difference model
mylm2 = lm(rsdiff ~ Cohabit, data = mydf)
summary(mylm2)


# ANCOVA of change model
mylm3 = lm(rsdiff ~ rsl + Cohabit, data = mydf)
summary (mylm3)


#### Data Generating Model: Full Equation

To illustrate the notion of two equations, let's rewrite the regression equation  

$$ score_{2i} = b_01_i + b_1score_{1i} + b_2cohabit_{i} + b_3(score_{1i} \times cohabit_{i}) + \epsilon_{i}$$

as two separate regression equations, one for fathers who graduated from highschool and one for fathers that did not. We can accomplish this by plugging in $0$ and $1$ into the regression equation and rearranging some of the terms. Doing so we get 

**Equation for Dyads That Did Cohabit**


$$ score_{2i} = (b_0 + b_2) +  (b_1 +  b_3)score_{1i} + \epsilon_{i}$$

**Equation for Dyads That Did Not Cohabit**

$$ score_{2i} = b_0 + b_1score_{1i} + \epsilon_{i}$$