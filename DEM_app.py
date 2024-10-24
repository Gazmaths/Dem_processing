pip --upgrade
import streamlit as st
import rasterio
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from rasterio.plot import show
from scipy.ndimage import gaussian_filter
from richdem import TerrainAttribute, LoadGDAL

# Upload files
st.title("Topographic Factor Extraction from DEM")
st.write("Upload a DEM (GeoTIFF) and a point shapefile for extracting topographic factors.")

dem_file = st.file_uploader("Upload DEM (GeoTIFF)", type=["tif"])
shapefile = st.file_uploader("Upload Shapefile (ZIP or SHP)", type=["zip", "shp"])

if dem_file and shapefile:
    # Read DEM
    with rasterio.open(dem_file) as dem:
        dem_data = dem.read(1)
        dem_affine = dem.transform
        dem_crs = dem.crs
        bounds = dem.bounds
        show(dem)

    # Display DEM properties
    st.write(f"DEM CRS: {dem_crs}")
    st.write(f"Bounds: {bounds}")
    
    # Compute Slope, Aspect, Curvature
    st.write("Processing DEM for topographic factors...")

    dem_rd = LoadGDAL(dem_file.name)  # Use RichDEM for topographic analysis
    
    # Calculate topographic factors
    slope = TerrainAttribute(dem_rd, attrib='slope_riserun')
    curvature = TerrainAttribute(dem_rd, attrib='curvature')
    twi = TerrainAttribute(dem_rd, attrib='twi')
    tpi = TerrainAttribute(dem_rd, attrib='tpi')

    # Visualize Slope
    st.subheader("Slope Map")
    fig, ax = plt.subplots()
    show(slope, ax=ax)
    st.pyplot(fig)

    # Load Shapefile and Extract Data
    if shapefile.name.endswith(".zip"):
        gdf = gpd.read_file(f"zip://{shapefile.name}")
    else:
        gdf = gpd.read_file(shapefile)

    st.write("Shapefile uploaded and read successfully!")

    # Reproject points to DEM CRS
    if gdf.crs != dem_crs:
        gdf = gdf.to_crs(dem_crs)

    # Extract topographic values for points
    st.write("Extracting topographic values for points...")
    points = gdf.geometry.apply(lambda x: (x.x, x.y))
    
    def extract_values(data, points, affine):
        indices = [~affine * point for point in points]
        return [data[int(row), int(col)] for row, col in indices]

    gdf["Slope"] = extract_values(slope, points, dem_affine)
    gdf["Curvature"] = extract_values(curvature, points, dem_affine)
    gdf["TWI"] = extract_values(twi, points, dem_affine)
    gdf["TPI"] = extract_values(tpi, points, dem_affine)

    st.write("Extraction Complete! Here's a preview:")
    st.write(gdf.head())

    # Allow Download of Results as CSV
    st.download_button("Download Extracted Data as CSV", gdf.to_csv(index=False), file_name="extracted_topographic_factors.csv")

else:
    st.write("Please upload both a DEM and a point shapefile.")
