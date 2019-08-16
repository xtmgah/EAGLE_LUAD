
################################################################################

rm(list = ls())
library(survival)
library(survey)

################################################################################

load("../data/data4Demo.RData")
#Input data: "data4Demo.RData" contains a data frame "dataPheno" and a list "distMatList"

#"dataPheno" is a subject by covariates table, it contains the variable needed in the survival analysis
#str(dataPheno)

#"distMatList" is a list of distance matrix, each component in this list is a n_i by n_i (n_i is the number of tumor samples of the i_th subject) symmetric matrix with all zeroes in diagonal
#These distance matrix can be calculated from any source of data (i.e. SCNA profiles, methylation profiles) by user's own definition
#str(distMatList)

#names of "distMatList" match the rownames of "dataPheno"
#table(rownames(dataPheno) == names(distMatList))

#column "nTumor" (number of tumor samples of that subject) in "dataPheno" matchs the size of distance matrix in "distMatList"
#table(dataPheno$nTumor == unlist(lapply(distMatList, nrow)))

################################################################################

#Calculate the APITH index from distMatList (average distance of each sample pair)
dataPheno$APITH <- unlist(lapply(distMatList, function(distMat)mean(distMat[lower.tri(distMat)])))
#hist(dataPheno$APITH, seq(0, 1, 0.05))

################################################################################

#Function to calculate weights for survival anlaysis using the inverse of the estimated variance of APITH index
#(Please refer to the supplementary notes for details of the algorithm in this function)

funWeight <- function(distMatList){
	tempFun <- function(distMat){
		k <- nrow(distMat)
		if(k <= 1)stop("Something is wrong!")
		distVec <- distMat[lower.tri(distMat)]
		k2 <- cbind(distVec, distVec)
		k3 <- NULL
		if(k >= 3){
			tempCombn <- combn(k, 3)
			for(i in 1:ncol(tempCombn)){
				tempIdx <- tempCombn[, i]
				tempMat <- cbind(c(distMat[tempIdx[1], tempIdx[2]], distMat[tempIdx[1], tempIdx[2]], distMat[tempIdx[1], tempIdx[3]]), c(distMat[tempIdx[1], tempIdx[3]], distMat[tempIdx[2], tempIdx[3]], distMat[tempIdx[2], tempIdx[3]]))
				k3 <- rbind(k3, tempMat)
			}
		}
		k4 <- NULL
		if(k >= 4){
			tempCombn <- combn(k, 4)
			for(i in 1:ncol(tempCombn)){
				tempIdx <- tempCombn[, i]
				tempMat <- cbind(c(distMat[tempIdx[1], tempIdx[2]], distMat[tempIdx[1], tempIdx[3]], distMat[tempIdx[1], tempIdx[4]]), c(distMat[tempIdx[3], tempIdx[4]], distMat[tempIdx[2], tempIdx[4]], distMat[tempIdx[2], tempIdx[3]]))
				k4 <- rbind(k4, tempMat)
			}
		}
		return(list(k2 = k2, k3 = k3, k4 = k4))
	}
	tempList <- lapply(distMatList, tempFun)
	k2 <- k3 <- k4 <- NULL
	for(i in 1:length(tempList)){
		k2 <- rbind(k2, tempList[[i]]$k2)
		k3 <- rbind(k3, tempList[[i]]$k3)
		k4 <- rbind(k4, tempList[[i]]$k4)
	}
	kk <- NULL
	for(i in 1:(length(distMatList) - 1)){
		for(j in (i + 1):length(distMatList)){
			tempDistMat1 <- distMatList[[i]]
			tempDistMat2 <- distMatList[[j]]
			tempDistVec1 <- tempDistMat1[lower.tri(tempDistMat1)]
			tempDistVec2 <- tempDistMat2[lower.tri(tempDistMat2)]
			tempMat <- cbind(rep(tempDistVec1, each = length(tempDistVec2)), rep(tempDistVec2, length(tempDistVec1)))
			kk <- rbind(kk, tempMat)
		}
	}
	a1 <- mean((k4[, 1] - k4[, 2])^2)
	a2 <- mean((k3[, 1] - k3[, 2])^2)
	a3 <- mean((kk[, 1] - kk[, 2])^2)
	b1 <- a1 / 2
	b2 <- (a1 - a2) / 2
	b3 <- (a3 - a1) / 2
	m <- unlist(lapply(distMatList, nrow))
	w <- 1 / (b3 + 2 / m /(m - 1) * b1 + 4 * (m - 2) / m / (m - 1) * b2)
	attr(w, "par") <- c(b1, b2, b3)
	return(w)
}

#Calculate the weights:
dataPheno$w <- funWeight(distMatList)
#str(dataPheno$w)

################################################################################

#Make "exp(coef)" as "odds ratio per 10% increase" of the APITH from SCNA
dataPheno$x <- dataPheno$APITH * 10

#Cox proportional hazard model without weight:
s1 <- summary(m1 <- coxph(SURVIVAL ~ STAGE + AGE_FIRST_DIAGNOSIS + SEX + x, data = dataPheno))
s2 <- summary(m2 <- coxph(METASTASIS ~ STAGE + AGE_FIRST_DIAGNOSIS + SEX + x, data = dataPheno))

#Cox proportional-hazard model without weight:with weight:
tempDesign <- svydesign(ids = ~ 1, data = dataPheno, weights = dataPheno$w)
s3 <- summary(m3 <- svycoxph(SURVIVAL ~ STAGE + AGE_FIRST_DIAGNOSIS + SEX + x, design = tempDesign, data = dataPheno))
s4 <- summary(m4 <- svycoxph(METASTASIS ~ STAGE + AGE_FIRST_DIAGNOSIS + SEX + x, design = tempDesign, data = dataPheno))

################################################################################

#Survival analysis of subjects with at least 3 tumor samples per subject
subDataPheno <- dataPheno[dataPheno$nTumor > 2, ]

#Cox proportional hazard model without weight:
s5 <- summary(m5 <- coxph(SURVIVAL ~ STAGE + AGE_FIRST_DIAGNOSIS + SEX + x, data = subDataPheno))
s6 <- summary(m6 <- coxph(METASTASIS ~ STAGE + AGE_FIRST_DIAGNOSIS + SEX + x, data = subDataPheno))

#Cox proportional-hazard model without weight:with weight:
tempDesign <- svydesign(ids = ~ 1, data = subDataPheno, weights = subDataPheno$w)
s7 <- summary(m7 <- svycoxph(SURVIVAL ~ STAGE + AGE_FIRST_DIAGNOSIS + SEX + x, design = tempDesign, data = subDataPheno))
s8 <- summary(m8 <- svycoxph(METASTASIS ~ STAGE + AGE_FIRST_DIAGNOSIS + SEX + x, design = tempDesign, data = subDataPheno))

################################################################################

# output tables

tempFun <- function(s)c(s$conf.int["x", c("exp(coef)", "lower .95", "upper .95")], P = s$coefficients["x", "Pr(>|z|)"])
tempWrite <- matrix(unlist(lapply(list(s1, s2, s3, s4, s5, s6, s7, s8), tempFun)), byrow = TRUE, ncol = 4)

dimnames(tempWrite) <- list(
	c("Survival_80subjects_no_weight", 
	"METASTASIS_80subjects_no_weight", 
	"Survival_80subjects_weight", 
	"METASTASIS_80subjects_weight", 
	"Survival_48subjects_no_weight", 
	"METASTASIS_48subjects_no_weight", 
	"Survival_48subjects_weight", 
	"METASTASIS_48subjects_weight"), 
	c("P", "OR", "Lower .95", "upper .95")
)
write.csv(tempWrite, file = "../table/table4Demo.csv", quote = FALSE)

################################################################################

# output figures

subDataPheno$x_cat3 <- cut(subDataPheno$APITH, c(0, 0.1, 0.3, 1))
#table(subDataPheno$x_cat3)

s <- summary(m <- coxph(SURVIVAL ~ STAGE + AGE_FIRST_DIAGNOSIS + SEX + strata(x_cat3), data = subDataPheno))

pdf("../figure/figure4Demo.pdf", width = 10, height = 10)
par(cex.lab = 2, cex.axis = 2, mar = c(5, 6, 2, 2))
plot(survfit(m), col = c("blue", "green", "red"), lwd = 2, xlab = "Survival Weeks", ylab = "Survival Rate")
legend("topright", c("Low ITH (16)", "Medium ITH (24)", "High ITH (8)"), lwd = 2, col = c("blue", "green", "red"), cex = 1.5)
dev.off()

################################################################################


