# GNSS‑InSAR 3D Deformation Estimation

> This project aims to achieve joint processing of multi‑source geodetic data using GNSS‑InSAR, supporting GNSS 3D ENU displacement, ascending/descending‑track D‑InSAR, ascending/descending‑track MAI‑InSAR observation fusion, LOS displacement 3D decomposition, multi‑source data variance component estimation, and other functions. All inputs are GeoTIFF raster data with geographic coordinates.

## Table of Contents

- [Introduction](#introduction)
- [Input Data Specification](#input-data-specification)
  - [File Overview](#file-overview)
  - [Data Preprocessing Requirements](#data-preprocessing-requirements)
- [Quick Start](#quick-start)
  - [Environment Dependencies](#environment-dependencies)
  - [Configure Input File Paths](#configure-input-file-paths)
  - [Run the Code](#run-the-code)
- [Dataset Access](#dataset-access)
- [Output Results](#output-results)
- [Contact](#contact)

## Introduction

This code processes GNSS 3‑D deformation, as well as ascending‑ and descending‑orbit D‑InSAR and MAI‑InSAR Line‑Of‑Sight (LOS) deformation observations, to carry out joint analysis of multi‑source geodetic measurements.

**Core Functions**

🟩 Read LOS displacement rasters from multi‑orbit D‑InSAR and MAI‑InSAR datasets;  
🟩 Automatically identify overlapping pixels among raster datasets;  
🟩 Automatically convert LOS displacements into 3‑D East‑North‑Up (ENU) components using radar incidence angles;  
🟩 Import accuracy rasters corresponding to each observation to serve as observation weights;  
🟩 Support GNSS rasters for the definition of 3‑D long‑wavelength errors;  
🟩 Export processed deformation products and quantitative evaluation metrics.

## Input Data Specification

> ⚠️ **Important**: All input files must be GeoTIFF with complete geospatial projection and spatial reference information. All rasters shall share an identical coordinate reference system (CRS) and pixel resolution. Invalid pixels must be assigned as `NaN`.

### File Overview

The core variables `Pline` and `Eline` in the code correspond to the filenames of ascending‑ and descending‑track D‑InSAR LOS‑displacement GeoTIFF files, respectively. The table below lists all example input files and their corresponding purposes:

| Example Filename | Orbit Type | Data Type | Description |
|---|---|---|---|
| `DA1239.tif` | Ascending | D‑InSAR LOS displacement | Ascending‑track D‑InSAR LOS raster (Corresponding variable: `Pline`) |
| `DASAR.tif` | Descending | D‑InSAR LOS displacement | Descending‑track D‑InSAR LOS raster (Corresponding variable: `Eline`) |
| `M1239.tif` | Ascending | MAI‑InSAR along‑track displacement | Ascending‑track MAI‑InSAR raster |
| `M0969.tif` | Descending | MAI‑InSAR along‑track displacement | Descending‑track MAI‑InSAR raster |
| `ILOS1239.tif` | Ascending | Radar incidence angle | Ascending‑track radar incidence angle (unit: degree), used for LOS displacement 3‑D decomposition |
| `ILOS0969.tif` | Descending | Radar incidence angle | Descending‑track radar incidence angle (unit: degree), used for LOS displacement 3‑D decomposition |
| `Precision1239.tif` | Ascending | D‑InSAR observation accuracy | Ascending‑track D‑InSAR observation accuracy raster, used as initial weight for observation |
| `Precision0969.tif` | Descending | D‑InSAR observation accuracy | Descending‑track D‑InSAR observation accuracy raster, used as initial weight for observation |
| `Precisionmai1239.tif` | Ascending | MAI‑InSAR observation accuracy | Ascending‑track MAI‑InSAR observation accuracy raster, used as initial weight for observation |
| `Precisionmai0969.tif` | Descending | MAI‑InSAR observation accuracy | Descending‑track MAI‑InSAR observation accuracy raster, used as initial weight for observation |
| `GPS1239e.tif` | — | GPS E‑displacement raster | GPS east‑displacement raster, used for defining 3‑D long‑wavelength errors |
| `GPS1239n.tif` | — | GPS N‑displacement raster | GPS north‑displacement raster, used for defining 3‑D long‑wavelength errors |
| `GPS1239u.tif` | — | GPS U‑displacement raster | GPS up‑displacement raster, used for defining 3‑D long‑wavelength errors |

### Data Preprocessing Requirements

🟩 All files shall be GeoTIFF with consistent geospatial reference;  
🟩 All rasters must have identical pixel resolution;  
🟩 Masked and invalid areas shall be assigned to `NaN` (NoData);  
🟩 Displacement units must be consistent across all datasets;  
🟩 Filenames of ascending‑/descending‑track MAI‑InSAR along‑track displacement files shall match those of ascending‑/descending‑orbit D‑InSAR displacement files;  
🟩 Filenames of ascending‑/descending‑track MAI‑InSAR along‑track precision files shall match those of ascending‑/descending‑orbit D‑InSAR displacement files;  
🟩 Filenames of ascending‑/descending‑track D‑InSAR LOS‑displacement precision files shall match those of ascending‑/descending‑orbit D‑InSAR displacement files;  
🟩 Filenames of ascending‑/descending‑track D‑InSAR radar incidence‑angle files shall match those of ascending‑/descending‑orbit D‑InSAR displacement files.

## Quick Start

### Environment Dependencies

🟩 MATLAB R2024b

### Configure Input File Paths

Modify the input‑related variables at the head of the script, and replace the example filenames with the **absolute paths** of your local files.

```matlab
% ====== User‑editable section: Input file configuration ======
Pline      = "1239";         % Filename prefix for ascending‑orbit D‑InSAR LOS displacement
> 💡 Note: It is recommended to remove any extra trailing space in the string variables to avoid runtime path errors.

### Run the Code

After completing the configuration, directly execute the main script `main.m`.

## Dataset Access

🟩 The raw measured datasets are large‑volume files and are **not included in this code repository**.
🟩 If you require the original test datasets, please contact the authors by leaving a message in GitHub Issues. You may also substitute the input files with your own processed InSAR and GNSS raster products following the format described above.

## Output Results

> 
> This is a template. Please revise it according to the actual outputs of your project. After execution finishes, outputs are saved to the current working directory by default.

Generated products include but are not limited to the following MATLAB variables:
% 3‑D deformation
X(1:size(prc,1)*3,1)

% Ascending‑track D‑InSAR long wavelength error
Aa*X(size(prc,1)*3+1:size(prc,1)*3+6,1)

% Descending‑track D‑InSAR long wavelength error
Ba*X(size(prc,1)*3+7:size(prc,1)*3+12,1)

% Ascending‑track MAI‑InSAR long wavelength error
Bma*Xma

% Descending‑track MAI‑InSAR long wavelength error
Bmd*Xmd

% Residual vector e for Ascending‑DInSAR, Descending‑DInSAR, Ascending‑MAI, Descending‑MAI, GNSS E‑N‑U
e = [Ya-[Aa Ba]*X(1:size(prc,1)*3+6,1); Yd-[Ad Bd]*Xd; Yma-[Ama Bma]*Xma; Ymd-[Amd Bmd]*Xmd; Yg-X(1:size(prc,1)*3,1)];
## Contact

🟩 Report bugs or submit suggestions: create GitHub Issues
🟩 Email: rakzfy@163.com
Eline      = "0969";        % Filename prefix for descending‑orbit D‑InSAR LOS displacement

