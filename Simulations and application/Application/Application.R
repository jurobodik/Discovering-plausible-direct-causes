####################    Application    ##########################
library("readxl")
library(countrycode)

GDP_dataset <- read.csv("GDP.csv")
head(GDP)

Education_dataset <- read.csv("Education.csv")
head(Education)

Mortality_and_Fertility <- read_excel("Mortality_and_Fertility.xlsx")
head(Mortality_and_Fertility)


data = list()
Countries = GDP_dataset$Country.Name
for (country in Countries) { if(!(country %in% Mortality_and_Fertility$`Region, subregion, country or area *`)){Countries = Countries[Countries!=country]}}


i=1
for (country in Countries) {
  GDP = t(  GDP_dataset[GDP_dataset$Country.Name == country,]  ); GDP = as.numeric( GDP[-1, ] )
  
  Education = t(  Education_dataset[Education_dataset$Country.Name == country,]  ); Education = as.numeric( Education[-1, ] )
  
  dataset_for_country = Mortality_and_Fertility[Mortality_and_Fertility$`Region, subregion, country or area *`==country,]
  dataset_for_country = dataset_for_country[dataset_for_country$Year>=1960,]
  dataset_for_country$GDP = GDP
  dataset_for_country$Education = Education
  
  data[[i]] = dataset_for_country
  i=i+1
}


continents_id=c()
continents = countrycode(sourcevar = Countries,
                         origin = "country.name",
                         destination = "continent")
continents[157] = "Africa"
for (i in 1:length(continents)) {
  if (continents[i] == "Americas") {continents_id[i] = 1 }
  if (continents[i] == "Africa") {continents_id[i] = 2 }
  if (continents[i] == "Oceania") {continents_id[i] = 3 }
  if (continents[i] == "Europe") {continents_id[i] = 4 }
  if (continents[i] == "Asia") {continents_id[i] = 5 }
}



X1 = c();
for (i in 1:length(Countries)) {
  X1 = c(X1, data[[i]]$GDP)
}

X2 = c();
for (i in 1:length(Countries)) {
  X2 = c(X2, data[[i]]$Education)
}

X3 = c();
for (i in 1:length(Countries)) {
  X3 = c(X3, data[[i]]$`Under-Five Mortality (deaths under age 5 per 1,000 live births)`)
}

X4 = c();
for (i in 1:length(Countries)) {
  X4 = c(X4, jitter(continents_id[i])) #jitter is here only because we use GAM. In the end, it fits the mean of each continent anyway
}

X5 = c();
for (i in 1:length(Countries)) {
  X5 = c(X5, data[[i]]$`Total Fertility Rate (live births per woman)...23`)
}

X = data.frame(as.numeric(X1), as.numeric(X2), as.numeric(X3), as.numeric(X4), as.numeric(X5) )

X = na.omit(X)

Y = X[,5]
X = X[,1:4]

set.seed(1);subset = sample(length(Y), 500) #You can try different set.seed, results will not change

Y = Y[subset]
X = X[subset,]
X1 = X[,1]; X2 = X[,2]; X3 = X[,3]; X4 = X[,4];

par(mfrow = c(3, 1))
plot(Y~X1, main = 'GDP')
plot(Y~X2, main = 'Education expenditure')
plot(Y~X3, main = 'Infant mortality rate')


#X1 = GDP
#X2 = Education
#X3 = Mortality
#X4 = Continent
#X5 = Fertility


IDS_and_score_based_estimation_of_F_parents(Y=Y, X=X, constraint_set_F = 'Location-scale') #change constraint_set_F = 'Gaussian', 'Additive', 'Pareto', 'Gamma' to see other results





