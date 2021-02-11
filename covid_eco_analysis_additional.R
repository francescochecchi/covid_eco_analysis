#..........................................................................................
###  ECOLOGICAL STUDY OF FACTORS ASSOCIATED WITH COVID-19 TRANSMISSIBILITY AND SEVERITY ###
#..........................................................................................

#..........................................................................................
## ----------- ADDITIONAL R CODE TO DOWNLOAD AND PREPARE CERTAIN DATASETS -------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Oct 2020)
                                          # (unless otherwise noted)

                                          # francesco.checchi@lshtm.ac.uk 


#..........................................................................................
### Preparatory steps
#..........................................................................................

  #...................................      
  ## Install or load required R packages
    
    # List of required packages
    x1 <- c("scales", "readxl", "data.table", "lubridate", "RColorBrewer", "gbm",
            "dismo", "conflicted", "gtools", "patchwork", "tidyverse", "countrycode", "MASS")
    
    # Install any packages not yet installed
    x2 <- x1 %in% row.names(installed.packages())
    if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
    
    # Load all packages    
    lapply(x1, library, character.only = TRUE)
    conflict_prefer("select", "dplyr")
    conflict_prefer("filter", "dplyr")
    conflict_prefer("area", "patchwork")


  #...................................      
  ## Starting steps

    # Clean up from previous code / runs
    rm(list=ls(all=TRUE) )
  
    # Set font
    windowsFonts(Arial=windowsFont("Arial"))

    # Set working directory to where this file is stored
    current_path = rstudioapi::getActiveDocumentContext()$path 
    setwd(dirname(current_path ))
    print( getwd() )
    
    # Initialise random numbers
    set.seed(123)
    
    # Colour-blind palette for graphing
    palette_cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    

#..........................................................................................
### Downloading or reading in required files
    
  #...................................
  ## Download R_t estimates by country over time developed by Imperial College MRC Outbreak Centre
      # Code written by Oliver Watson and colleagues (Imperial College)

    # Function to download estimates
    latest_rt <- function() {

      # parse for latest date
      url <- "https://github.com/mrc-ide/global-lmic-reports/tree/master/data"
      html <- xml2::read_html(url)
      links <- rvest::html_nodes(html, ".js-navigation-open.link-gray-dark")
      text <- rvest::html_text(links, "title")
      latest <- which.max(as.Date(head(substr(text, 1, 10), -1)))
      d_link <- paste0("https://github.com/mrc-ide/global-lmic-reports/raw/master/data/",text[latest])

      # download and read in data file
      tf <- tempfile()
      download.file(d_link, tf)
      tf2 <- tempfile()
      extract <- unzip(tf, exdir = tf2)
      dat <- as.data.frame(data.table::fread(extract))

      # subset to Rt values
      dat <- dat[dat$compartment == "Rt" & dat$scenario == "Maintain Status Quo",]

      return(dat)
    }

    # Download and write file
    out <- latest_rt()
    write.csv(out, "rt_imperial.csv", row.names = FALSE)


  #...................................
  ## Download Google mobility change estimates by country over time, while averaging mobility estimates across different environments
      # Code written by Oliver Watson and colleagues (Imperial College)

    # Function to download while also averaging across different mobility environments
    download_url <- function(url) {
      tryCatch({
        tf <- tempfile()
        code <- download.file(url, tf, mode = "wb")
        if (code != 0) {
          stop("Error downloading file")
        }
      },
      error = function(e) {
        stop(sprintf("Error downloading file '%s': %s, please check %s",
                     url, e$message))
      })
      return(tf)
    }


    # Load Google Mobility Data
    date <- as.Date(date)
    goog_tf <- download_url("https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv")
    goog <- read.csv(goog_tf, stringsAsFactors = FALSE)
    goog$iso3c <- countrycode::countrycode(goog$country_region, "country.name", "iso3c")

      # N.B. The date format changed recently in Google and may change again. Look out for errors related to this
      # the subregions change over time so catch for what should be included if this breaks.
      # we just want to filter to country level measures

    # Select and average data of interest
    mob <-  goog %>%
      filter(sub_region_1 == "" & sub_region_2 == "" & metro_area == "") %>%
      mutate(overall = 1/4 * retail_and_recreation_percent_change_from_baseline +
               1/4 * grocery_and_pharmacy_percent_change_from_baseline +
               1/4 * transit_stations_percent_change_from_baseline +
               1/4 * workplaces_percent_change_from_baseline) %>%
      mutate(date = as.Date(date, format="%Y-%m-%d")) %>%
      select(country_region, iso3c, date, overall)

    # Write file
    write.csv(mob, "mobility.csv", row.names = FALSE)

