# TITLE: PISCO_data_Prep
# AUTHOR: Andres Pinos-Sanchez

# This is script is to put together the dataset to be used as part of the final 
# project of Bayesian stats class (OSU), and to estimate parameter values for the
# urchin section of the ORTM project for Andres' thesis project

# O.G. data: https://data.piscoweb.org/metacatui/view/doi%3A10.6085%2FAA%2Fpisco_recruitment.1477.1

#------------------------------- TOP LEVEL STUFF ------------------------------

# Load libraries
library(ggplot2)
library(dplyr)
library(readr)

# Set working directory (Why not)
setwd("C:/Users/pinosa/OneDrive - Oregon State University/Ed_OSU_Thesis_GradSchool/OR_Trophic_Model/Data/PISCO_Recruitment_Data")

# - - - - - - - - - - - - - - - End of section - - - - - - - - - - - - - - - -



#------------------- PISCO (all) DATA EXPLORATION & CLEANING ------------------

# Input data from PISCO, data name corresponds to collection year (not deployment)
PISCO_data_1989 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.124.1.csv")  |> mutate(year = 1989, .before = 1)
PISCO_data_1990 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.126.1.csv")  |> mutate(year = 1990, .before = 1)
PISCO_data_1991 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.128.1.csv")  |> mutate(year = 1991, .before = 1)
PISCO_data_1992 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.130.1.data") |> mutate(year = 1992, .before = 1)
PISCO_data_1993 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.132.1.data") |> mutate(year = 1993, .before = 1)
PISCO_data_1994 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.134.1.csv")  |> mutate(year = 1994, .before = 1)
PISCO_data_1995 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.136.1.data") |> mutate(year = 1995, .before = 1)
PISCO_data_1996 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.138.1.data") |> mutate(year = 1996, .before = 1)
# 1997 skipped
PISCO_data_1998 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.140.1.data") |> mutate(year = 1998, .before = 1)
PISCO_data_1999 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.142.1.data") |> mutate(year = 1999, .before = 1)
PISCO_data_2000 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.144.1.csv")  |> mutate(year = 2000, .before = 1)
PISCO_data_2001 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.146.1.data") |> mutate(year = 2001, .before = 1)
PISCO_data_2002 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.148.1.csv")  |> mutate(year = 2002, .before = 1)
PISCO_data_2003 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.150.1.data") |> mutate(year = 2003, .before = 1)
PISCO_data_2004 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.152.1.csv")  |> mutate(year = 2004, .before = 1)
PISCO_data_2005 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.154.1.csv")  |> mutate(year = 2005, .before = 1)
PISCO_data_2006 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.156.1.data") |> mutate(year = 2006, .before = 1)
PISCO_data_2007 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.1364.1.data")|> mutate(year = 2007, .before = 1)
PISCO_data_2008 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.1366.1.csv") |> mutate(year = 2008, .before = 1)
PISCO_data_2009 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.1438.1.csv") |> mutate(year = 2009, .before = 1)
PISCO_data_2010 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.1446.1.csv") |> mutate(year = 2010, .before = 1)
PISCO_data_2011 <- read_csv("Datasets/doi_10.6085_AA_pisco_recruitment.1464.1.txt") |> mutate(year = 2011, .before = 1)

# # Print and explore
# print()
# summary()

# - - - - - - - - - - - - - - - End of section - - - - - - - - - - - - - - - -



#----------------------------- BIND & COMBINE DATA ----------------------------

# Combine all datasets
PISCO_all_years <- bind_rows(PISCO_data_1989, PISCO_data_1990, PISCO_data_1991, 
                             PISCO_data_1992, PISCO_data_1993, PISCO_data_1994, 
                             PISCO_data_1995, PISCO_data_1996, PISCO_data_1998,
                             PISCO_data_1999, PISCO_data_2000, PISCO_data_2001,
                             PISCO_data_2002, PISCO_data_2003, PISCO_data_2004, 
                             PISCO_data_2005, PISCO_data_2006, PISCO_data_2007, 
                             PISCO_data_2008, PISCO_data_2009, PISCO_data_2010, 
                             PISCO_data_2011)

# # Print and explore
# print()
# summary()
 
# Save table
write.csv(PISCO_all_years, "PISCO_all_years.csv", row.names = FALSE)

# - - - - - - - - - - - - - - - End of section - - - - - - - - - - - - - - - -
