library(nlme)     # CRAN v3.1-163 
library(survival) # CRAN v3.5-7 
library(JM)       # CRAN v1.5-2 

# Load in data
DF <- read.csv("wearablesensorstudy.csv")

# Longitudinal submodel
lmefit <- lme(VMI ~ TIME + LP + LP:TIME,
  random = ~ TIME + I(TIME^2) | ID, data = DF,
  control = lmeControl(opt = "optim")
)

DS <- DF[!rev(duplicated(rev(DF$ID))), ]

# Time-to-event submodel
coxfit <- coxph(Surv(STIME, STATUS) ~ LP,
  data = DS, x = TRUE
)

# Fit joint model
JMfit <- jointModel(lmeObject = lmefit, survObject = coxfit, timeVar = "TIME")
summary(JMfit)
