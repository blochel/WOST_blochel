
# climate data 2 ----------------------------------------------------------

# https://www.erikkusch.com/courses/krigr/quickstart/ 


library('KrigR')
library('terra')
library('tidyverse')


sessionInfo()



API_User <- "alexanderblochel@ufl.edu"
API_Key <- "9601e245-c02f-415a-914b-1c30f11b3a5b"


Dir.Base <- setwd('C:/Users/alexanderblochel/OneDrive - University of Florida/Desktop/WOST_nesting_initiative/data/climate_data/') # identifying the current directory
Dir.Data <- file.path(Dir.Base, "data") # folder path for data
Dir.Covariates <- file.path(Dir.Base, "Covariates") # folder path for covariates
Dir.Exports <- file.path(Dir.Base, "Exports") # folder path for exports
## create directories, if they don't exist yet
Dirs <- sapply(
  c(Dir.Data, Dir.Covariates, Dir.Exports),
  function(x) if (!dir.exists(x)) dir.create(x)
)
rm(Dirs) # we don't need to keep the response to directory creation


Extent_ext <- ext(c(-86.308594, -70.090820, 23.805450, 32.175612))
Extent_ext

#dates 
start_date <- '2014-01-01 00:00'
end_date <- '2023-12-31 24:00'

## Note that I have already downloaded the QuickStart raw data and CDownloadS() is simply loading this from the disk for me here. Your console output while CDownloadS() is being executed will look differently.
QuickStart_Raw <- CDownloadS(
  ## Variable and Data Product
  Variable = "2m_temperature", # this is air temperature
  DataSet = "reanalysis-era5-land", # data product from which we want to download
  ## Time-Window
  DateStart = start_date, # date at which time window opens
  DateStop = end_date, # date at which time window terminates
  TZone = "GMT", # European Central Time to align with our study region
  ## Temporal Aggregation
  TResolution = "day", # we want daily aggregates
  TStep = 1, # we want aggregates of 1 day each
  ## Spatial Limiting
  Extent = Extent_ext, # our rectangular bounding box
  ## File Storing
  Dir = Dir.Data, # where to store the data
  FileName = 'temperature_raw_FL2', # what to call the resulting file
  ## API User Credentials
  API_User = API_User,
  API_Key = API_Key
)


Plot.SpatRast(QuickStart_Raw[[110]])


metags(QuickStart_Raw) %>% 
  filter(name == 'Citation')


## Note that I have already downloaded the global GMTED2010 data with this function prior, your output will show the download itself as well
Covs_ls <- CovariateSetup(
  Training = QuickStart_Raw,
  Target = .02,
  Dir = Dir.Covariates,
  Keep_Global = TRUE
)
#The CovariateSetup() function can also be used to prepare raster data you already have at hand for use in subsequent Kriging. 

Covs_ls



QuickStart_Krig_FL <- Kriging(
  Data = QuickStart_Raw, # data we want to krig as a raster object
  Covariates_training = Covs_ls[[1]], # training covariate as a raster object
  Covariates_target = Covs_ls[[2]], # target covariate as a raster object
  Equation = "GMTED2010", # the covariate(s) we want to use
  nmax = 40, # degree of localisation
  Cores = parallel::detectCores()-1, # we want to krig using as many cores possible to speed this process up
  FileName = "QuickStart_Krig_FL2", # the file name for our full kriging output
  Dir = Dir.Exports, # which directory to save our final input in
  Compression = 9, #compress to speed up 
  Keep_Temporary = FALSE, 
  verbose = TRUE
)

Plot.SpatRast(QuickStart_Krig_FL[[110]])
plot(QuickStart_Krig_FL)

terra::plot(QuickStart_Krig_FL$Prediction[1])


QuickStart_Krig_FL$Prediction %>% 
  as.data.frame()

#C:/Users/alexanderblochel/OneDrive - University of Florida/Desktop/WOST_nesting_initiative/scripts/Climate_data/Exports/QuickStart_Krig_FL2_Kriged.nc

sessionInfo()

