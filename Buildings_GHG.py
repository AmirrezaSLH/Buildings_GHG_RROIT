# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 16:26:16 2024

@author: asalehi
"""

import pandas as pd
import geopandas as gpd

import os
import numpy as np

import matplotlib.pyplot as plt

# Initialize Model Data
def init_emissions_factors(file_path='Database/Default/Emissions_Factors.csv'):
    """
    Initialize emissions factors from a CSV file.
    
    This function reads a CSV file containing emissions factors for various sources across different years.
    It organizes these factors into a dictionary for easy access.

    Parameters:
    - file_path (str): The file path to the CSV file containing emissions factors.

    Returns:
    - dict: A dictionary where each key is an emissions source and its value is another dictionary 
            mapping years to emissions factors for that source.
    """
    # Read the CSV file into a DataFrame
    emissions_factors_df = pd.read_csv(file_path)
    
    # Define the column name for the year
    year_column = 'YEAR'
    # Identify all emission source columns (excluding the year column)
    emissions_sources = set(emissions_factors_df.columns) - {year_column}
    
    # Initialize a dictionary to store emissions factors by source and year
    emissions_factors_dict = {}
    
    # Populate the dictionary with emissions factors for each source
    for emissions_source in emissions_sources:
        emissions_factor_source_dict = {}
        
        for idx, row in emissions_factors_df.iterrows():
            year = int(row[year_column])
            emissions_factor = row[emissions_source]
            emissions_factor_source_dict[year] = emissions_factor
        
        emissions_factors_dict[emissions_source] = emissions_factor_source_dict
    
    return emissions_factors_dict

def init_btype_info(file_path='Database/Default/BTYPE_info.csv'):
    """
    Initialize building type information from a CSV file.
    
    This function reads a CSV file containing energy use intensity (EUI) information for various building types.
    It organizes this information into a dictionary for easy access and returns a list of building types.

    Parameters:
    - file_path (str): The file path to the CSV file containing building type information.

    Returns:
    - tuple:
        - list: A list of building types code(BTYPE_CODE).
        - dict: A dictionary where each key is a building type (BTYPE_CODE) and its value is another dictionary
                mapping EUI categories to their values for that building type.
    """
    # Read the CSV file into a DataFrame
    btype_info_df = pd.read_csv(file_path)
    
    # Extract a list of building types (BTYPE_CODE) from the DataFrame
    bt_list = btype_info_df["BTYPE_CODE"].tolist()
    
    # Define the column name for the building type
    building_type_column = 'BTYPE_CODE'
    building_type_fullname_column = 'BTYPE'
    # Identify all EUI columns (excluding the building type column)
    eui_columns = set(btype_info_df.columns) - {building_type_column} - {building_type_fullname_column}
    
    # Initialize a dictionary to store EUI information by building type
    all_bt_eui_dict = {}
    
    # Populate the dictionary with EUI information for each building type
    for idx, row in btype_info_df.iterrows():
        btype = row[building_type_column]
        bt_eui_dict = {et: row[et] for et in eui_columns}
        all_bt_eui_dict[btype] = bt_eui_dict
    
    return bt_list, all_bt_eui_dict


def init_degree_days(file_path='Database/default/Degree_Days.csv'):
    """
    Initialize Heating Degree Days (HDD) and Cooling Degree Days (CDD) from a CSV file.
    
    This function reads a CSV file containing degree day information (HDD and CDD) for various years 
    and organizes them into two dictionaries: one for HDD and one for CDD.
    
    Parameters:
    - file_path (str): The file path to the CSV file containing degree day information. 
                       Defaults to 'Database/default/Degree_Days.csv'.
    
    Returns:
    - tuple: A tuple containing two dictionaries:
        - HDD_dict (dict): A dictionary where the keys are years and the values are HDD values.
        - CDD_dict (dict): A dictionary where the keys are years and the values are CDD values.
    """
    
    # Read the CSV file containing degree day data into a DataFrame
    degree_days_df = pd.read_csv(file_path)
    
    # Initialize empty dictionaries to store HDD and CDD data by year
    HDD_dict = {}
    CDD_dict = {}
    
    # Iterate through each row in the DataFrame to populate the HDD and CDD dictionaries
    for idx, row in degree_days_df.iterrows():
        year = row['YEAR']  # Extract the year
        HDD = row['HDD']    # Extract the Heating Degree Days (HDD) for that year
        CDD = row['CDD']    # Extract the Cooling Degree Days (CDD) for that year
        
        # Populate the dictionaries with year as the key and HDD/CDD as the values
        HDD_dict[year] = HDD
        CDD_dict[year] = CDD
    
    # Return the two dictionaries: one for HDD and one for CDD
    return HDD_dict, CDD_dict

#USER

def init_bt_mapping(filepath='Database/User/BT_mapping.csv'):
    """
    Reads a Building Type (BT) mapping from a CSV file and creates a dictionary mapping from user BT codes to model BTs.
    
    This function loads a CSV file that maps user-defined building type codes (USER_BT_CODE) to model building types (MODEL_BT_CODE).
    It converts this mapping into a dictionary for easy access and use in further processing or mapping tasks.

    Parameters:
    - filepath (str, optional): The file path to the CSV file containing the BT mappings. Defaults to a predefined path.

    Returns:
    - dict: A dictionary where each key is a user-defined building type code (USER_BT_CODE) and its value is the corresponding
            model building type (MODEL_BT_CODE).
    """
    # Load the CSV file into a DataFrame
    bt_mapping_df = pd.read_csv(filepath)
    
    # Directly convert the DataFrame to a dictionary
    bt_mapping_dict = bt_mapping_df.set_index('USER_BT_CODE')['MODEL_BT_CODE'].to_dict()
    
    return bt_mapping_dict

#RROIT Data
    
def init_buildings_shapefile(file_path='Database/RROIT/Buildings/Buildings_shapefile.shp'):
    """
    Read a buildings shapefile, ensure the CRS is EPSG:26917, calculate the footprint area, and filter columns.
    
    This function loads a shapefile into a GeoDataFrame, checks if the CRS is EPSG:26917 or sets it to this value if not,
    calculates the footprint area of each geometry, and retains only specified columns in the final GeoDataFrame.
    
    Parameters:
    - file_path (str): The file path to the shapefile containing building geometries.

    Returns:
    - GeoDataFrame: A GeoDataFrame with the CRS set to EPSG:26917, an additional 'footprint' column for the area,
                    and filtered to include only 'Build_Type', 'Land_Use', 'footprint', and 'geometry' columns.
    """
    # Load the shapefile into a GeoDataFrame
    buildings_shapefile = gpd.read_file(file_path)
    #remove empty rows
    buildings_shapefile = buildings_shapefile[buildings_shapefile['geometry'].notna()]
    
    # Check and adjust the CRS to EPSG:26917
    if buildings_shapefile.crs and buildings_shapefile.crs.to_string() == 'EPSG:26917':
        print("The CRS of the GeoDataFrame is EPSG:26917.")
    elif buildings_shapefile.crs:
        print(f"Original CRS is set: {buildings_shapefile.crs}")
        buildings_shapefile = buildings_shapefile.to_crs(epsg=26917)  # Convert the CRS to EPSG:26917
        print(f"Converted to CRS: {buildings_shapefile.crs}")
    else:
        print("CRS is not set. Setting to EPSG:26917.")
        buildings_shapefile.crs = "EPSG:26917"  # Assuming the intention to set to EPSG:26917 if CRS is not present
    
    #Set ID
    buildings_shapefile['ID'] = buildings_shapefile.index
    
    # Calculate the footprint area of each building and add it as a new column
    buildings_shapefile['footprint'] = buildings_shapefile['geometry'].area
    
    #Find centroid of each building
    buildings_shapefile['centroid'] = buildings_shapefile.geometry.centroid
    
    # Filter the GeoDataFrame to keep only specified columns
    columns_to_keep = ['ID', 'Build_Type', 'Land_Use', 'footprint', 'centroid', 'geometry']
    
    buildings_shapefile = buildings_shapefile[columns_to_keep]
    
    return buildings_shapefile

def init_landuse_shapefile(file_path='Database/RROIT/Landuse/Landuse_shapefile.shp'):
    """
    Read a landuse shapefile, ensure the CRS is EPSG:26917, calculate the footprint area, and filter columns.
    
    This function loads a landuse shapefile into a GeoDataFrame, checks if the CRS is EPSG:26917 or sets it to this value if not,
    calculates the footprint area of each geometry, and retains only specified columns in the final GeoDataFrame.
    
    Parameters:
    - file_path (str): The file path to the shapefile containing landuse geometries.

    Returns:
    - GeoDataFrame: A GeoDataFrame with the CRS set to EPSG:26917, an additional 'footprint' column for the area,
                    and filtered to include only 'FID', 'TYPE', 'footprint', and 'geometry' columns.
    """
    # Load the shapefile into a GeoDataFrame
    landuse_shapefile = gpd.read_file(file_path)
    #remove empty rows
    landuse_shapefile = landuse_shapefile[landuse_shapefile['geometry'].notna()]
    
    # Check and adjust the CRS to EPSG:26917
    if landuse_shapefile.crs and landuse_shapefile.crs.to_string() == 'EPSG:26917':
        print("The CRS of the GeoDataFrame is EPSG:26917.")
    elif landuse_shapefile.crs:
        print(f"Original CRS is set: {landuse_shapefile.crs}")
        landuse_shapefile = landuse_shapefile.to_crs(epsg=26917)  # Convert the CRS to EPSG:26917
        print(f"Converted to CRS: {landuse_shapefile.crs}")
    else:
        print("CRS is not set. Setting to EPSG:26917.")
        landuse_shapefile.crs = "EPSG:26917"  # Assuming the intention to set to EPSG:26917 if CRS is not present
    
    #Set ID
    landuse_shapefile['ID'] = landuse_shapefile.index
    
    # Calculate the footprint area of each landuse and add it as a new column
    landuse_shapefile['footprint'] = landuse_shapefile['geometry'].area
    
    # Filter the GeoDataFrame to keep only specified columns
    columns_to_keep = ['ID', 'FID', 'TYPE', 'footprint', 'geometry']
    landuse_shapefile = landuse_shapefile[columns_to_keep]
    
    return landuse_shapefile

def init_land_building_mapping(bt_mapping_dict, file_path='Database/RROIT/Landuse/BuildingsLanduse_CV_SWMP.csv'):
    """
    Initialize a mapping between land use types and building types based on percentages from a CSV file.
    
    This function reads a CSV file containing mappings between building codes and land use types, along with a
    percentage that represents the distribution of each building type within a land use type. It uses a provided
    dictionary to map building codes to model building types codes (MODEL_BT_CODE), aggregates the percentages by land use
    and MODEL_BT_CODE, pivots the table for a comprehensive view, and finally converts it into a dictionary of dictionaries
    for easy access.
    
    Parameters:
    - bt_mapping_dict (dict): A dictionary mapping from building codes to model building types.
    - file_path (str, optional): The file path to the CSV file containing the mappings. Defaults to a predefined path.

    Returns:
    - dict: A dictionary where each key is a land use type and its value is another dictionary mapping MODEL_BT_CODE to
            the aggregated percentage sum of that building type for the land use.
    """
    # Read the CSV file into a DataFrame
    land_building_mapping_df = pd.read_csv(file_path)
    
    # Map 'buildingCode' to the model building type (MODEL_BT_CODE) using the provided dictionary
    land_building_mapping_df['MODEL_BT_CODE'] = land_building_mapping_df['buildingCode'].map(bt_mapping_dict)
    
    # Group by 'landuseType' and 'MODEL_BT_CODE', and sum 'percent' for each group
    aggregated_df = land_building_mapping_df.groupby(['landuseType', 'MODEL_BT_CODE'])['percent'].sum().reset_index()
    
    # Pivot the aggregated DataFrame to have 'MODEL_BT_CODE' as columns with their sums as values
    # Fill missing values with 0 to ensure every MODEL_BT_CODE has a value for each landuseType
    pivot_df = aggregated_df.pivot(index='landuseType', columns='MODEL_BT_CODE', values='percent').fillna(0)
    
    # Convert the pivoted DataFrame to a dictionary of dictionaries for land use to building type ratios
    landuse_building_ratio_dict = pivot_df.to_dict(orient='index')
    
    return landuse_building_ratio_dict

#Spatial Manupilations

def landuse_contain_buildings(buildings_shapefile, landuse_shapefile, landuse_building_ratio_dict):
    """
    Filters out landuse entries from the landuse_shapefile that do not contain buildings or 
    have landuse types not listed in the landuse_building_ratio_dict. It iterates through each 
    landuse entry, checks its type against the provided dictionary, and removes it if not found 
    or if the type is None.

    Parameters:
    - buildings_shapefile (GeoDataFrame): Not used in the current implementation, but intended for future use.
    - landuse_shapefile (GeoDataFrame): A GeoDataFrame containing land use types.
    - landuse_building_ratio_dict (dict): A dictionary with land use types as keys.

    Returns:
    - GeoDataFrame: The filtered landuse_shapefile with only the relevant land use entries.
    """
    
    # Extract the list of valid land use types from the dictionary keys
    landuse_list = landuse_building_ratio_dict.keys()
    
    # Iterate over each row in the landuse_shapefile
    for idx, row in landuse_shapefile.iterrows():
        # Retrieve the land use type for the current row
        landuse_type = row['TYPE']
        
        # Check if the landuse type is either not in the valid list or is None
        if landuse_type not in landuse_list:
            # Drop the row if its landuse type is not in the list
            landuse_shapefile = landuse_shapefile.drop(idx)
            print(f"landuse index {idx} type is {landuse_type} and contains no building, dropped.")
        elif landuse_type == None:
            # Drop the row if its landuse type is None
            landuse_shapefile = landuse_shapefile.drop(idx)
            print(f"landuse index {idx} type is None, insufficient information for further analysis, dropped.")
        
    return landuse_shapefile

def buildings_landuse_mapping(buildings_shapefile, landuse_shapefile):
    """
    Maps buildings to land use types based on the spatial location of building centroids within land use areas.
    
    This function ensures each building in the buildings_shapefile has a centroid calculated. It then performs
    a spatial join between these centroids and the landuse_shapefile to determine the land use category for each building.
    The resulting land use information (ID and type) is added back to the original buildings shapefile.
    
    Parameters:
    - buildings_shapefile (GeoDataFrame): A GeoDataFrame containing geometries of buildings.
    - landuse_shapefile (GeoDataFrame): A GeoDataFrame containing land use areas with their geometries and types.
    
    Returns:
    - GeoDataFrame: The updated buildings_shapefile with 'landuse_ID' and 'landuse_type' columns added.
    """

    # Ensure the 'centroid' column exists in buildings_shapefile
    if 'centroid' not in buildings_shapefile.columns:
        buildings_shapefile['centroid'] = buildings_shapefile.geometry.centroid
    
    # Convert centroids to a new GeoDataFrame for spatial join
    centroids = gpd.GeoDataFrame(geometry=buildings_shapefile['centroid'], crs=buildings_shapefile.crs)
    
    # Perform spatial join between centroids and selected columns of landuse_shapefile
    # This operation matches each building centroid with its corresponding land use area
    joined = gpd.sjoin(centroids, landuse_shapefile[['geometry', 'ID', 'TYPE']], how='left', op='within')
    
    # Remove duplicates to ensure each building is mapped to only one land use type
    joined = joined.drop_duplicates(subset='geometry', keep='first')
    
    # Add the resulting land use information back to the original buildings shapefile
    buildings_shapefile['landuse_ID'] = joined['ID']
    buildings_shapefile['landuse_type'] = joined['TYPE']
    
    return buildings_shapefile

# GHG Calculations

def EUI_weather_norm(all_bt_eui_dict, year, HDD_dict, CDD_dict, avg_HDD, avg_CDD):
    """
    Normalize Energy Use Intensity (EUI) values based on weather data, specifically Heating Degree Days (HDD)
    and Cooling Degree Days (CDD), for a given year.

    Parameters:
    - all_bt_eui_dict (dict): A dictionary where keys are building types, and values are dictionaries mapping 
                              fuel types (e.g., 'NG', 'ELECTRICITY') to EUI values.
    - year (int): The year for which the normalization is to be applied.
    - HDD_dict (dict): A dictionary mapping years to Heating Degree Days (HDD).
    - CDD_dict (dict): A dictionary mapping years to Cooling Degree Days (CDD).
    - avg_HDD (float): The average Heating Degree Days for normalization.
    - avg_CDD (float): The average Cooling Degree Days for normalization.

    Returns:
    - dict: A normalized EUI dictionary where the values are adjusted based on weather data (HDD/CDD).
    """
    
    # Initialize an empty dictionary to store the normalized EUI values
    all_bt_eui_normalized_dict = {}
    
    # Iterate through each building type in the EUI dictionary
    for bt in all_bt_eui_dict.keys():
        all_bt_eui_normalized_dict[bt] = {}  # Initialize an empty dictionary for each building type
        
        # Iterate through each fuel type for the current building type
        for fuel_type in all_bt_eui_dict[bt].keys():
            
            # Normalize 'Electricity' usage based on CDD
            if fuel_type == 'ELECTRICITY':
                all_bt_eui_normalized_dict[bt][fuel_type] = (all_bt_eui_dict[bt][fuel_type] / avg_CDD) * CDD_dict[year]
                
            # Normalize 'Natural Gas (NG)' usage based on HDD
            elif fuel_type == 'NG':
                all_bt_eui_normalized_dict[bt][fuel_type] = (all_bt_eui_dict[bt][fuel_type] / avg_HDD) * HDD_dict[year]
            
            #Normalize other fuels that are usually used for heating
            else:
                all_bt_eui_normalized_dict[bt][fuel_type] = (all_bt_eui_dict[bt][fuel_type] / avg_HDD) * HDD_dict[year]
                
            
    # Return the dictionary with normalized EUI values
    return all_bt_eui_normalized_dict



def bt_emissions_intensity_calculator(all_bt_eui_dict, emissions_factors_dict, year, HDD_dict, CDD_dict, avg_HDD, avg_CDD):
    """
    Calculates the emissions intensity for different building types based on their Energy Use Intensity (EUI),
    weather-normalized EUI values for a specific year, and emissions factors for various fuel types.

    Parameters:
    - all_bt_eui_dict (dict): A dictionary where keys are building types and values are dictionaries mapping
                              fuel types to their EUI values.
    - emissions_factors_dict (dict): A dictionary where keys are fuel types and values are dictionaries mapping
                                     years to emissions factors for that fuel type.
    - year (int): The year for which to calculate emissions intensity.

    Returns:
    - dict: A dictionary where keys are building types and values are their calculated emissions intensity for the given year.
    """
    
    # Normalize EUI values for the specified year
    all_bt_eui_normalized_dict = EUI_weather_norm(all_bt_eui_dict, year, HDD_dict, CDD_dict, avg_HDD, avg_CDD)

    # Initialize a dictionary to store the calculated emissions intensity for each building type
    bt_emissions_intensity_dict = {}
    
    # Iterate over each building type in the normalized EUI dictionary
    for bt in all_bt_eui_normalized_dict.keys():
        bt_eui_dict = all_bt_eui_normalized_dict[bt]
        bt_emissions_intensity = 0  # Initialize emissions intensity for the building type
        
        # Iterate over each fuel type in the building type's EUI dictionary
        for fuel_type in bt_eui_dict.keys():
            fuel_type_eui = bt_eui_dict[fuel_type]  # EUI for the fuel type
            fuel_type_emissions_factor = emissions_factors_dict[fuel_type][year]  # Emissions factor for the fuel type
            
            # Calculate emissions intensity for the fuel type and add it to the building type's total
            fuel_type_emissions_intensity = fuel_type_emissions_factor * fuel_type_eui
            bt_emissions_intensity += fuel_type_emissions_intensity
        
        # Store the calculated emissions intensity for the building type
        bt_emissions_intensity_dict[bt] = bt_emissions_intensity
    
    return bt_emissions_intensity_dict

def spatial_GHG_calculator_building_level(buildings_shapefile, land_building_mapping, bt_emissions_intensity_dict, bt_mapping):
    """
    Calculates the greenhouse gas (GHG) emissions for each building in the shapefile.
    The calculation is based on the building's footprint, its land use type, and the emissions intensity
    of its building type. The results are stored in a new "GHG" column within the buildings shapefile.

    Parameters:
    - buildings_shapefile (GeoDataFrame): The shapefile containing buildings data.
    - land_building_mapping (dict): A mapping between land use types and building types, including their probabilities.
    - bt_emissions_intensity_dict (dict): A dictionary mapping building types to their emissions intensity.
    - bt_mapping (dict): A dictionary mapping building types in the shapefile to those used in the emissions intensity dictionary.

    Returns:
    - GeoDataFrame: The input shapefile with an added "GHG" column representing calculated GHG emissions for each building.
    """

    BGHG_list = []  # List to hold the calculated GHG values for each building
    
    # Iterate over each building in the shapefile
    for index, row in buildings_shapefile.iterrows():
        BGHG = 0  # Initialize the GHG variable for the building
        
        footprint = row['footprint']  # Building footprint area
        
        landuse = row['landuse_type']  # Land use type of the building
        landuse_list = land_building_mapping.keys()  # List of all possible land use types
        
        building_type = row['Build_Type']  # Original building type from the shapefile
        # Map the building type to a known type if possible; otherwise, set to None
        building_type_mapped = bt_mapping.get(building_type, None)
        building_type_list = bt_emissions_intensity_dict.keys()  # List of all building types in the emissions intensity dict
        
        #Checks if the dataset has specific GHG assignments
        if 'GHG' in row.index:
            if pd.notna(row['GHG']):
                BGHG = row['GHG'] 
            
        # Check if the building has a known  building type
        elif building_type_mapped in building_type_list:
            emission_intensity = bt_emissions_intensity_dict[building_type_mapped]
            BGHG = footprint * emission_intensity  # Calculate GHG based on footprint and emissions intensity
   
        # If the building type is not known check if it has a known landuse
        elif landuse in landuse_list:
            # Calculate GHG for buildings with indirect mappings through land use
            for bt in building_type_list:
                emission_intensity = bt_emissions_intensity_dict[bt]
                probability = (land_building_mapping[landuse].get(bt, 0))/100  # Default probability is 0 if not explicitly set
                print(bt_emissions_intensity_dict[bt])
                BGHG += probability * emission_intensity * footprint  # Aggregate GHG based on probability and emissions intensity
                
        else:
            # Set GHG to None if no suitable mapping is found
            BGHG = None

        BGHG_list.append(BGHG)  # Append the calculated GHG value to the list
    
    # Add the GHG calculations as a new column in the shapefile
    #kg CO2e
    buildings_shapefile["GHG"] = BGHG_list
    # Remove any buildings that do not have a GHG calculation
    buildings_shapefile = buildings_shapefile.dropna(subset=["GHG"])
    
    #Calculate GHG intensity
    #kg CO2e / m2
    buildings_shapefile['GHG_Intensity'] = buildings_shapefile["GHG"] / buildings_shapefile["footprint"]
    
    #Calculate GHG in Tonnes
    buildings_shapefile["GHG_T"] = buildings_shapefile["GHG"] / 1000
    return buildings_shapefile


def spatial_GHG_calculator_landuse_level(buildings_shapefile, landuse_shapefile):
    """
    Aggregates GHG emissions at the landuse level by summing the emissions of individual buildings
    within each landuse category. The aggregated GHG data is then merged into the landuse shapefile.

    Parameters:
    - buildings_shapefile (GeoDataFrame): A GeoDataFrame containing buildings data, including their GHG emissions.
    - landuse_shapefile (GeoDataFrame): A GeoDataFrame containing landuse data.
    - Land_BT_ratio_dic (dict): Not used in the current implementation, but intended for future enhancements
                                regarding the ratio of building types within each land use category.
    - BTGHG_dic (dict): Not used in the current implementation, but intended for future enhancements
                        regarding GHG emissions data by building type.

    Returns:
    - GeoDataFrame: The updated landuse shapefile with an additional 'GHG' column representing the aggregated
                    GHG emissions for each land use category.
    """

    # Aggregate GHG emissions for each landuse ID in the buildings shapefile
    ghg_aggregated = buildings_shapefile.groupby('landuse_ID')['GHG'].sum()
    
    # Convert the Series object to a DataFrame for merging
    ghg_aggregated_df = ghg_aggregated.reset_index(name='GHG')
    
    # Merge the aggregated GHG data with the landuse shapefile based on landuse ID
    #kg CO2e
    landuse_shapefile = landuse_shapefile.merge(ghg_aggregated_df, left_on='ID', right_on='landuse_ID', how='left')
    
    #kg CO2e / m2
    landuse_shapefile['GHG_Intensity'] = landuse_shapefile["GHG"] / landuse_shapefile["footprint"]
    
    # GHG in Tonnes
    landuse_shapefile["GHG_T"] = landuse_shapefile["GHG"]/1000
    
    return landuse_shapefile

def initialize_building_GHG_layers():
    """
    Initialize and prepare all necessary data for calculating greenhouse gas (GHG) emissions for buildings.
    
    This function loads various data sources, including emissions factors, building type information, shapefiles
    for buildings and land use, and mappings between building types and land uses. It filters and processes
    the data to prepare for GHG calculations.
    
    Returns:
    - bt_info (dict): Detailed information for each building type.
    - emissions_factors (dict): Emissions factors for various fuel types.
    - buildings_shapefile_mapped (GeoDataFrame): A GeoDataFrame of buildings mapped to their corresponding land use.
    - land_building_ratio (dict): A mapping of land use types to building types and their ratios.
    - bt_mapping (dict): A mapping of actual building types to model building types.
    - landuse_shapefile_filtered (GeoDataFrame): A filtered GeoDataFrame containing only relevant land use data.
    """

    # Initialize emissions factors from a predefined source, which provides fuel-specific emissions factors.
    emissions_factors = init_emissions_factors()
    
    # Initialize building type information, which includes a list of building types and detailed EUI data for each type.
    bt_list, bt_info = init_btype_info()
    
    # Initialize shapefiles for both buildings and land uses from predefined sources.
    buildings_shapefile = init_buildings_shapefile()
    landuse_shapefile = init_landuse_shapefile()

    # Initialize a mapping from actual building types (as recorded in the shapefile) to model building types.
    bt_mapping = init_bt_mapping()

    # Initialize a mapping from land use types to building types, including their respective ratios.
    land_building_ratio = init_land_building_mapping(bt_mapping)
    
    # Filter the land use shapefile to include only relevant land uses that contain buildings.
    landuse_shapefile_filtered = landuse_contain_buildings(buildings_shapefile, landuse_shapefile, land_building_ratio)
    
    # Map buildings to their corresponding land use types, adding this information to the buildings shapefile.
    buildings_shapefile_mapped = buildings_landuse_mapping(buildings_shapefile, landuse_shapefile_filtered)
    
    # Return all the necessary initialized and processed data for further GHG calculations.
    return bt_info, emissions_factors, buildings_shapefile_mapped, land_building_ratio, bt_mapping, landuse_shapefile_filtered
 
def run_buildings_GHG(bt_info, emissions_factors, buildings_shapefile_mapped, land_building_ratio, bt_mapping, landuse_shapefile_filtered, year):
    """
    Run the process to calculate and plot greenhouse gas (GHG) emissions for buildings and land uses for a given year.

    This function calculates emissions intensity based on provided building type information, emissions factors, and weather data.
    It then calculates GHG emissions at both the building and land use levels.

    Parameters:
    - bt_info (dict): Detailed building type information.
    - emissions_factors (dict): Emissions factors for different fuel types.
    - buildings_shapefile_mapped (GeoDataFrame): Shapefile containing building geometries, footprints, and mapped land use data.
    - land_building_ratio (dict): A mapping between land use types and building types, including their ratios.
    - bt_mapping (dict): A mapping from real-world building types to model building types.
    - landuse_shapefile_filtered (GeoDataFrame): Filtered shapefile of land use areas containing buildings.
    - year (int): The year for which GHG emissions are to be calculated.

    Returns:
    - tuple: 
        - GeoDataFrame with GHG emissions at the building level.
        - GeoDataFrame with GHG emissions aggregated at the land use level.
    """

    # Initialize Degree Days (HDD and CDD) from a predefined source
    HDD_dict, CDD_dict = init_degree_days()
    
    # Average Degree Days for normalization (can be adjusted based on your dataset)
    avg_HDD, avg_CDD = 3265.64, 486.88  # Placeholder values, to be adjusted with actual data These are the average values between 2021- 2050 from Climate Atlas Canada
    
    # Normalize the building energy use intensity (EUI) based on weather data (Degree Days)
    bt_emissions_intensity = bt_emissions_intensity_calculator(bt_info, emissions_factors, year, HDD_dict, CDD_dict, avg_HDD, avg_CDD)

    temp_building_shapefile = buildings_shapefile_mapped.copy()
    # Calculate GHG emissions for each building based on its footprint, type, and land use
    buildings_GHG = spatial_GHG_calculator_building_level(temp_building_shapefile, land_building_ratio, bt_emissions_intensity, bt_mapping)
    
    temp_landuse_shapefile = landuse_shapefile_filtered.copy()
    # Aggregate GHG emissions at the land use level by summing the building-level GHG data within each land use area
    landuse_GHG = spatial_GHG_calculator_landuse_level(buildings_GHG, temp_landuse_shapefile)
    
    # Return the calculated GHG emissions at both the building and land use levels
    return buildings_GHG, landuse_GHG

#%%

#The following functions are for saving the outputs and creating plots.

def create_output_file(building_ghg_gdf, landuse_ghg_gdf, building_output_name="bghg", landuse_output_name='lghg', path="output", output_type="ESRI Shapefile"):
    """
    Saves GeoDataFrames to specified file formats after processing.

    Parameters:
    - building_ghg_gdf (GeoDataFrame): GeoDataFrame containing building GHG data.
    - landuse_ghg_gdf (GeoDataFrame): GeoDataFrame containing landuse GHG data.
    - building_output_name (str): Base name for the building output file.
    - landuse_output_name (str): Base name for the landuse output file.
    - path (str): Directory path where the files will be saved.
    - output_type (str): Type of the output file ('ESRI Shapefile', 'GPKG', 'GeoJSON').
    """

    # Create the directory if it doesn't exist
    if not os.path.exists(path):
        os.makedirs(path)
    
    # Prepare building GHG data for export
    building_ghg_gdf_copy = building_ghg_gdf.copy()
    building_ghg_gdf_copy.drop('centroid', axis=1, inplace=True)
    building_ghg_gdf_copy.rename(columns={'GHG_Intensity': 'GHGIntense', 'landuse_type': 'LU_Type'}, inplace=True)
    
    # Prepare landuse GHG data for export
    landuse_ghg_gdf_copy = landuse_ghg_gdf.copy()
    landuse_ghg_gdf_copy.rename(columns={'GHG_Intensity': 'GHGIntense'}, inplace=True)
    
    # Define file paths
    fpath_building = os.path.join(path, building_output_name)
    fpath_landuse = os.path.join(path, landuse_output_name)
    extension = {'ESRI Shapefile': '.shp', 'GPKG': '.gpkg', 'GeoJSON': '.geojson'}
    
    # Check if the output type is supported and prepare the final path
    if output_type in extension:
        fpath_building_final = fpath_building + extension[output_type]
        fpath_landuse_final = fpath_landuse + extension[output_type]

        # Save files in the specified format
        building_ghg_gdf_copy.to_file(fpath_building_final, driver=output_type)
        landuse_ghg_gdf_copy.to_file(fpath_landuse_final, driver=output_type)
        
        print(f'Files saved successfully in {output_type} format.')
    else:
        print(f'Unsupported file format: {output_type}')

    return 0

def plot_heatmap(gdf, column_name, plot_title, legend_title, cmap='viridis', figsize=(10, 6)):
    """
    Plot a GeoDataFrame with color mapping based on a specified column.

    Parameters:
    gdf : GeoDataFrame
        GeoDataFrame containing the spatial data and values to plot.
    column_name : str
        Name of the column to be used for the color mapping.
    plot_title : str
        Title of the plot.
    legend_title : str
        Title of the legend, representing the units or the data type.
    cmap : str, optional
        Colormap for the plot (default is 'viridis').
    figsize : tuple, optional
        Size of the figure (default is (10, 6)).
    """
    # Create figure and axes
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    
    # Plot the GeoDataFrame
    gdf_plot = gdf.plot(column=column_name, ax=ax, legend=True, cmap=cmap, legend_kwds={'label': legend_title})
    
    # Add the title to the plot
    ax.set_title(plot_title, fontsize=15)
    
    # Remove the frame of the plot
    ax.set_frame_on(False)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
    # Move legend to the right side of the plot
    
    '''
    cbar = fig.colorbar(gdf_plot.collections[0], ax=ax)
    cbar.set_label(legend_title)
    '''
    '''
    leg = ax.get_legend()
    if leg:  # Check if legend is there (in case column_name is not valid, it won't be)
        leg.set_bbox_to_anchor((1, 0.5))
        leg.set_title(legend_title)
        '''
    # Display the plot
    plt.show()
    
    return 0

#%%

"""
The following code is an example of how to use the module.

It demonstrates the process on Cooksville data.
"""

# Initialize the model with all necessary data layers.
# Run this line once. After the first run, the model will be initialized, and this step doesn't need to be repeated.
bt_info, emissions_factors, buildings_shapefile_mapped, land_building_ratio, bt_mapping, landuse_shapefile_filtered = initialize_building_GHG_layers()


# Calculate GHG emissions for a specific year (e.g., 2020).
# This step can be run multiple times for different years to get year-specific results.
bghg, lghg = run_buildings_GHG(bt_info, emissions_factors, buildings_shapefile_mapped, land_building_ratio, bt_mapping, landuse_shapefile_filtered, 2020)

# You can repeat the line above with different years if needed, e.g., for 2021:

#%%    

"""
Demonstration: Plotting GHG Emission Maps and Printing Aggregate Results
"""

# Plot heatmap for building-level total GHG emissions
plot_heatmap(bghg, 'GHG_T', 'Cooksville Buildings GHG Emissions', 'tonnes $CO_{2}e$', cmap='jet', figsize=(10, 6))

# Plot heatmap for building-level GHG intensity (per square meter)
plot_heatmap(bghg, 'GHG_Intensity', 'Cooksville Buildings GHG Emissions Intensities', 'Kg $CO_{2}e$/$M^{2}$', cmap='jet', figsize=(10, 6))

# Plot heatmap for land use-level total GHG emissions
plot_heatmap(lghg, 'GHG_T', 'Cooksville Landuse GHG Emissions', 'tonnes $CO_{2}e$', cmap='jet', figsize=(10, 6))

# Plot heatmap for land use-level GHG intensity (per square meter)
plot_heatmap(lghg, 'GHG_Intensity', 'Cooksville Landuse GHG Emission Intensities', 'Kg $CO_{2}e$/$m^{2}$', cmap='jet', figsize=(10, 6))


# Print Aggregate Emission Estimates
print("Total Building GHG Emissions (kg CO2e):", bghg['GHG'].sum())  # Total kg CO2e
print("Total Building GHG Emissions (Tonnes CO2e):", bghg['GHG_T'].sum())  # Total tonnes CO2e
print("Total Building GHG Emissions (Kilo Tonnes CO2e):", bghg['GHG_T'].sum() / 1000) 
