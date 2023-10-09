# P1_Yield.R
# About: This program will open the county yield and climate data the preproccess
#       the data. Once the data has been preprocesses a very simple regression 
#       will be applied. 
#     
# Inputs:  CornYield, climte data
# Outputs: 
#
# T. A. Schillerberg
#               Sep. 2023
#      Updated: Oct. 2023

# Office Computer
# setwd("C:/Users/tas0053/OneDrive - Auburn University/Research/DataScienceExp/Script/")
# fileloc1 <- 'C:/Users/tas0053/OneDrive - Auburn University/Research/DataScienceExp/Data/'
# Mac
setwd("~/Library/CloudStorage/OneDrive-AuburnUniversity/Research/DataScienceExp/Script")
fileloc1 <- '~/Library/CloudStorage/OneDrive-AuburnUniversity/Research/DataScienceExp/Data/'

# Libraries ####################################################################
library(tidyverse)
library(MASS)
# library(caret)
# library(zoo)

# Part I Variables To Change ###################################################
loc1 <- c('NASS/','NOAA/')
yr <- 1950:2022

# Part II Functions ############################################################
naSum <- function(dat){
  sum(is.na(dat))
}
detrend <- function(dat, yr){
  MASS::rlm(dat ~ yr, method="M")$residuals
}
lmYield <- function(dat){
  datF <- tibble('DatY' = t(dat[1:length(yr)]), 
                'DatX1' = t(dat[(length(yr)+1):(2*length(yr))]),
                'DatX2' = t(dat[(2*length(yr)+1): (3*length(yr))]),
                'DatX3' = t(dat[(3*length(yr)+1): (4*length(yr))]),
                'DatX4' = t(dat[(4*length(yr)+1): (5*length(yr))]))
  colnames(datF) <- c('DatY', 'DatX1', 'DatX2', 'DatX3', 'DatX4')
  lm1 <- lm(formula = DatY ~ DatX1 + DatX2 + DatX3 + DatX4, data = datF,
            na.action = na.omit)
  result <- lm1$coefficients %>% unlist()
  return(result)
}
# Part III Exploring Yield #####################################################
# . 3.1 Opening and combining yield --------------------------------------------
datCorn <- rbind(read_csv(paste0(fileloc1,loc1[1],'CornYield_19502022_ALKS','.csv'),
                      col_names = TRUE, col_types = cols('c')),
             read_csv(paste0(fileloc1,loc1[1],'CornYield_19502022_KYNY','.csv'),
                      col_names = TRUE, col_types = cols('c')),
             read_csv(paste0(fileloc1,loc1[1],'CornYield_19502022_NCVA','.csv'),
                      col_names = TRUE, col_types = cols('c')),
             read_csv(paste0(fileloc1,loc1[1],'CornYield_19502022_WAWY','.csv'),
                      col_names = TRUE, col_types = cols('c')))
# We read the data in as a character, because we want to maintain the '0' for 
# the state and county fips codes. As this will be needed later on when we are 
# combining the climate data. Note in the summary not all columns were read in 
# as a character.
summary(datCorn)
head(datCorn)
dim(datCorn)

# . 3.2 Yield Pre-processing ---------------------------------------------------
# Removing rows of "Other Counties"
datCorn <- datCorn %>%
  filter(!County == 'OTHER COUNTIES' ) %>%
  filter(!County == 'OTHER (COMBINED) COUNTIES')
# Removing unwanted columns
datCorn <- datCorn %>%
  subset(select = -c(Program, Period, `Week Ending`, `Geo Level`, `Ag District`, 
                     `Ag District Code`, `Zip Code`, Region, watershed_code, 
                     Watershed, Commodity, `Data Item`, Domain, 
                     `Domain Category`, `CV (%)`))
dim(datCorn)

# Changing from a long
datCorn <- spread(datCorn, key = Year, value = Value, fill = NA)
dim(datCorn)
summary(datCorn)

# Plot of number of available data by year
y1 <- apply(datCorn[,5:ncol(datCorn)], MARGIN = 2, FUN = naSum)
dat <- tibble(yr, y1)
plot(dat)

# Remove the counties with more than 30% of their data missing
x1 <- apply(datCorn[,5:ncol(datCorn)], MARGIN = 1, FUN = naSum)
x1 <- x1/length(1950:2022)
datCorn[x1 > 0.3,] <- NA
datCorn <- datCorn[!is.na(datCorn$State),]
dim(datCorn)

# Plot of number of available data by year
y1 <- apply(datCorn[,5:ncol(datCorn)], MARGIN = 2, FUN = naSum)
dat <- tibble(yr, y1)
plot(dat)

# # Filling missing values
# head <- datCorn[,1:4] 
# dat <- caret::preProcess(datCorn[,5:ncol(datCorn)], method = 'knnImpute',
#                          k = 10)
# datCorn <- t(datCorn) %>%
#   cbind(yr)
# colnames(datCorn) <- c(paste0(head$`State ANSI`,head$`County ANSI`), 'Year')

# Remove the rows with missing years
datCorn <- na.omit(datCorn)
dim(datCorn)
head <- datCorn[,1:4] 
datCorn <- t(datCorn[5:ncol(datCorn)]) %>%
  cbind(yr)
# datCorn <- as_tibble(datCorn)
colnames(datCorn) <- c(paste0(head$`State ANSI`,head$`County ANSI`), 'Year')

# De-trend the data using linear with M estimator
ggplot(as_tibble(datCorn), aes(x=yr, y = `01049`)) +
  geom_line()
datCorn <- apply(datCorn[,datCorn[1:(ncol(datCorn)-1)]], MARGIN = 2, 
                 FUN = detrend, yr = yr) %>%
  as_tibble() %>%
  cbind(yr)
colnames(datCorn) <- c(paste0(head$`State ANSI`,head$`County ANSI`), 'Year')
ggplot(as_tibble(datCorn), aes(x=yr, y = `01049`)) +
  geom_line()

# Part IV Exploring Climate Data ###############################################
# . 4.1 Opening and combining yield --------------------------------------------
datP <- read_fwf(paste0(fileloc1,loc1[2],'ncei.noaa.gov_pub_data_cirs_climdiv',
                        '_climdiv-pcpncy-v1.0.0-20230907',
                        '.txt'),
                 col_types = cols('c'))
datD <- read_fwf(paste0(fileloc1,loc1[2],'ncei.noaa.gov_pub_data_cirs_climdiv_',
                        'climdiv-pdsicy-v1.0.0-20230907',
                        '.txt'),
                 col_types = cols('c'))
datTx <- read_fwf(paste0(fileloc1,loc1[2],'ncei.noaa.gov_pub_data_cirs_climdiv_',
                         'climdiv-tmaxcy-v1.0.0-20230907',
                        '.txt'),
                 col_types = cols('c'))
datTn <- read_fwf(paste0(fileloc1,loc1[2],'ncei.noaa.gov_pub_data_cirs_climdiv_',
                         'climdiv-tmincy-v1.0.0-20230907',
                        '.txt'),
                 col_types = cols('c'))

summary(datP); summary(datD); summary(datTx); summary(datTn)
head(datP)
dim(datP); dim(datD); dim(datTx); dim(datTn)
# From the readme file of the climate data the state code is in position 1-2, 
# county FIPS 3-5, Element code 6-7, and year 8-11
# . 4.2 Climate Pre-processing -------------------------------------------------
# Correct data type in datTx
datTx$X7 <- as.numeric(datTx$X7)
datTx$X8 <- as.numeric(datTx$X8)
datTx$X9 <- as.numeric(datTx$X9)
summary(datTx)

# taking care of NAs according to documentation
summary(datP)
datP[datP == -9.99] <- NA
summary(datP)
datD[datD == -99.99] <- NA
datTx[datTx == -99.99] <- NA
datTn[datTn == -99.99] <- NA

# Separating out the State and County Code
datStateP <- apply(datP[,1], MARGIN = 1, FUN= substr, 1,5)
yearP <- apply(datP[,1], MARGIN = 1, FUN= substr, 8,11)
datStateD <- apply(datD[,1], MARGIN = 1, FUN= substr, 1,5)
yearD <- apply(datD[,1], MARGIN = 1, FUN= substr, 8,11)
datStateTx <- apply(datTx[,1], MARGIN = 1, FUN= substr, 1,5)
yearTx <- apply(datTx[,1], MARGIN = 1, FUN= substr, 8,11)
datStateTn <- apply(datTn[,1], MARGIN = 1, FUN= substr, 1,5)
yearTn <- apply(datTn[,1], MARGIN = 1, FUN= substr, 8,11)

# Create new data frame of only the growing season month average
grow <- apply(datP[,5:10], MARGIN = 1, FUN = sum, na.rm=TRUE)
datP <- tibble(StateCo = datStateP, Year = yearP, Precip = grow)
grow <- apply(datD[,5:10], MARGIN = 1, FUN = mean, na.rm=TRUE)
datD <- tibble(StateCo = datStateD, Year = yearD, PDSI = grow)
grow <- apply(datTx[,5:10], MARGIN = 1, FUN = mean, na.rm=TRUE)
datTx <- tibble(StateCo = datStateTx, Year = yearTx, Tmax = grow)
grow <- apply(datTn[,5:10], MARGIN = 1, FUN = mean, na.rm=TRUE)
datTn <- tibble(StateCo = datStateTn, Year = yearTn, Tmin = grow)

# Create a Wide data frame
dim(datP)
datP <- spread(datP, key = Year, value = Precip)
datD <- spread(datD, key = Year, value = PDSI)
datTx <- spread(datTx, key = Year, value = Tmax)
datTn <- spread(datTn, key = Year, value = Tmin)

# Select only the rows that we need
dim(datP)
m <- which(datP$StateCo %in% colnames(datCorn))
datP <- datP[m,]
m <- which(datD$StateCo %in% colnames(datCorn))
datD <- datD[m,]
m <- which(datTx$StateCo %in% colnames(datCorn))
datTx <- datTx[m,]
m <- which(datTn$StateCo %in% colnames(datCorn))
datTn <- datTn[m,]
dim(datP); dim(datD); dim(datTx); dim(datTn)
identical(datP$StateCo, datD$StateCo) 
identical(datP$StateCo,datTx$StateCo)
identical(datP$StateCo, datTn$StateCo)

# Select only the years that we need
datP <- datP %>%
  subset(select = c(StateCo, `1950`:`2022`))
datD <- datD %>%
  subset(select = c(StateCo, `1950`:`2022`))
datTx <- datTx %>%
  subset(select = c(StateCo, `1950`:`2022`))
datTn <- datTn %>%
  subset(select = c(StateCo, `1950`:`2022`))

# Part V Regression ############################################################
# . 5.1 Formatting -------------------------------------------------------------
# Select the yield counties that match the climate data
dim(datCorn)
m <- which(colnames(datCorn) %in% datP$StateCo)
datCorn <- datCorn[,m] %>%
  t()
dim(datCorn)

# Make sure that the order is the same
order <- rownames(datCorn)
datP <- datP %>%
  arrange(factor(StateCo, levels = order))
datD <- datD %>%
  arrange(factor(StateCo, levels = order))
datTx <- datTx %>%
  arrange(factor(StateCo, levels = order))
datTn <- datTn %>%
  arrange(factor(StateCo, levels = order))


# Fitting the regression
row <- 35
datS <- rbind(datCorn[row,], datP[row,2:ncol(datP)], datD[row,2:ncol(datD)], 
             datTx[row,2:ncol(datTx)], datTn[row,2:ncol(datTn)]) %>% 
  t() %>%
  as_tibble()
dim(datS)
colnames(datS) <- c('DatY', 'DatX1', 'DatX2', 'DatX3', 'DatX4')
lm(formula = DatY ~ DatX1 + DatX2 + DatX3 + DatX4, data = datS,
   na.action = na.omit) %>%
  summary()

dat <- cbind(datCorn, datP[,2:ncol(datP)], datD[,2:ncol(datD)], 
             datTx[,2:ncol(datTx)], datTn[,2:ncol(datTn)])
result <- matrix(data = 0, nrow = dim(dat)[1], ncol = 5)
for (i in 1:dim(datCorn)[1]){
  x <- lmYield(dat[i,]) %>% unlist()
  result[i,] <- x
}

result <- as_tibble(result)
colnames(result) <- c('Intercept','Precip.','PDSI','Tmax','Tmin')
result <- cbind('StateCo' = order, result)



# END ##########################################################################