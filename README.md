# EGM722_Assessment
Repository for EGM722 module assessment

# WIND FARM SITE SUITABILITY ANALYSIS IN NORTHERN IRELAND

## 1 INTRODUCTION
The purpose of this repository is to assist in wind farm site suitability analysis in Northern Ireland (NI). Several data files are provided to assess the suitability of sites based on site characteristics and impacts to humans and the environment in terms of their proximity to the sites. Several proposed sites are passed to the scripts/notebook producing a locator map, two CSVs (one with proximity analysis results and one with site characteristics) and two interactive map HTMLs (one of which is optional) in the output_files folder. The locator map is created to display the proposed sites in the context of NI as a whole. The proximity analysis and site characteristics CSV outputs allow for detailed comparisons between the sites and the proximity analysis can be rerun at different buffer distances. The interactive map(s) allow for the investigation of transport networks to be used during wind farm construction but can also be used to further visualise the proximity and site characteristic results as the CSV data is displayed in popups on the map(s).

## 2 SETUP AND INSTALLATION
### 2.1 REQUIRED INSTALLATIONS
Firstly, the user should have git installed on their computer. Installation instructions can be found [here](https://git-scm.com/downloads). Secondly, the user should have conda installed. Anaconda Navigator is the GUI for conda and installation instructions can be found [here](https://docs.anaconda.com/anaconda/install/). 

The user should check system requirements and install the appropriate version for their device’s operating system. 
### 2.2 RECOMMENDED INSTALLATIONS	
It is recommended to install GitHub Desktop (the user will need to login using their GitHub account). Installation instructions can be found [here](https://desktop.github.com/). GitHub Desktop is a GUI for git and GitHub, among other things, it allows the user to save changes locally before committing them to their remote repository.

The user should use an IDE to run the code within this repository. The code was created using PyCharm and it is therefore recommended to use this IDE if possible. PyCharm can be downloaded and installed using instructions found [here](https://www.jetbrains.com/pycharm/download/), be sure to download the *community* edition as this is free. 
### 2.3 FORKING THE REPOSITORY
Once all installations are complete the user can go ahead and copy (fork) this repository to their GitHub account by clicking ‘fork’ in the top right-hand corner of this window making a note of the new repository URL.
### 2.4 CLONING THE REPOSITORY
The user should now clone their copied repository to their device. Two suggested ways of doing this are:
1.	Open GitHub Desktop, select ‘File’ > ‘Clone Repository’. Select the URL tab and enter the user’s repository URL (which will likely be https://<GitHub username>/EGM722_Assessment). Save to a local folder on the user’s device (and make a note of this folder path). Click ‘Clone’. When asked how you plan to use this fork in the next window select ‘For my own purposes’. 
2.	Open Git Bash (from the start menu). Navigate to the local folder the user wishes to clone their repository to and execute ‘git clone <repository URL>.git’. The user should see some messages regarding downloading/unpacking files and the repository should then be set up. 
### 2.5 CREATING/SETTING UP THE ENVIRONMENT
Open Anaconda Navigator, select the ‘Environments’ tab on the left-hand side. Select ‘Import’, this option can be found at the bottom left of the page. The user should navigate to the folder path on their local drive where they cloned their repository and select the ‘environments.yml’ file provided. This file contains the environment name, channels to install packages from and a list of main dependencies (with some versions specified for compatibility). <br>
  
**Name:** egm722_assessment<br>
**Channels:** - conda-forge - defaults<br>
**Dependencies:** - python=3.9 - geopandas - cartopy>=0.21 - notebook=6.5.3 - rasterio - pyepsg - folium - numpy - rasterstats>=0.18<br>
  
Ensure an appropriate environment name is being used and select ‘Import’.  After some time, all the required dependencies should have been imported. 
### 2.6 PYCHARM SETUP (ONLY IF USING PYCHARM)
Open PyCharm and select ‘New Project’. For ‘Location’ – the user should navigate to the folder path on their local drive where they cloned their repository. Under ‘Python Interpreter’ select ‘Previously configured interpreter’, ‘Add Interpreter’ and ‘Add Local Interpreter’. On the window that opens the user should select ‘Conda Environment’ and ‘Use existing environment’ and navigate to the python interpreter that is part of their conda environment and their conda executable (program) (this is relative to where they’ve installed Anaconda).  Click ‘Create’ and if a pop up appears select ‘Create from Existing Sources’. 

If the user already uses PyCharm select ‘File’ in the top-left then ‘New Project’ and follow the remaining instruction as above. 
### 2.7 INPUT DATA
Data files are provided within the repository data_files folder covering NI. Files containing proposed wind farm site locations (Site_Locations) are also provided for illustrative purposes. The user is directed to the ‘Data_sources_and_descriptions’ file within the data_files folder for more information on the provided input data including the uniquely identifying naming column labels. 

The user should upload their proposed site locations into the data_files folder (ensure the files are named Site_Locations and the uniquely identifying naming column is called ‘Name’). 

The user can update any/all input files to improved/updated data. The user can also use input data for regions other than NI adding data of the same type as contained within the repository but for their study area and/or adding new data types not covered. All new input data should obey the following:
-	Use the same file names as the files already present within the repository
-	Be saved in the data_files folder
-	Vector data should contain a naming column that uniquely identifies the data (where possible matching those of the data already supplied*)
-	Data must be spatially referenced (no naive geometries)
-	Place information CSV files must have a Name, Type, County, LGD and Population column<br>

**If using different unique naming column labels to those in the provided files the user will need to update the naming column labels within all code scripts.*
### 2.8 RUNNING THE SCRIPTS
The user should open pycharm and jupyter notebook from the anaconda navigator home page with their correct environment activated. They should then navigate to the repository scripts/notebook. The user should update the scripts as per section 2.8.1 below and then run them in the following order:
  1. Locator_map.py
  2. Proximity_analysis.py
  3. Site_characteristics.py
  4. Interactive_map.ipynb
#### 2.8.1 UPDATES REQUIRED BEFORE RUNNING THE SCRIPTS
Before running the code, the user should update the length of the colour lists to match the number of sites within their newly uploaded Site_Locations. The colour lists are found in the Locator_map.py (line 293) and Interactive_map.ipynb (9th code cell down). Additional colours to be used can be found [here](https://matplotlib.org/stable/gallery/color/named_colors.html).

If the user increases or decreases the number of unique values within the Counties name column or Place_information type column, much like when adding Site_Locations, the user must update the length of the display lists within the Locator_map.py script (lines 239, 264, 268-9, 272 and 274). 

Also, if using this repository for another region other than NI, the user must update the CRS values at the start of the Locator_map.py, Proximity_analysis.py and Site_characteristics.py scripts after the input data files are loaded. Please bear in mind that the scripts are designed for a spatial reference system that has map units in metres. 

Line 198 in the Locator_map.py script clips the lakes data to the extent of the study preventing them from appearing outside the study outline on the map. If the user updates their input data, they can either remove this part of the script or reproduce it for other data they wish to clip to the study extent outline.

The Proximity analysis.py script is set to run with a buffer of 5km however the user can update this before running the analysis or to repeat the analysis for multiple buffer distances. The user should define their desired buffer distance in code lines 351-352.

## 3 TROUBLESHOOTING
The user may experience errors when running the scripts. If these occur when running a function, the user is advised to show the function docstring by calling help(), the docstring will explain the required inputs and data types. The user is directed to the ‘Troubleshooting’ file in this repository for more information on some errors they may encounter and how to fix them.
