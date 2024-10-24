import streamlit as st
import rasterio
import numpy as np
import geopandas as gpd
from rasterio.plot import show
import matplotlib.pyplot as plt

# File Upload Section
st.title("Topographic Analysis App")
dem_file = st.file_uploader("Upload DEM (GeoTIFF format)", type=["tif", "tiff"])
shapefile = st.file_uploader("Upload Point Shapefile", type=["shp", "dbf", "shx", "prj"])

def calculate_tpi(dem, window_size=3):
    # Calculate TPI using a moving window (mean of neighbors)
    kernel = np.ones((window_size, window_size))
    kernel[window_size // 2, window_size // 2] = 0
    local_mean = np.pad(dem, pad_width=window_size//2, mode='reflect')
    local_mean = (local_mean[1:, 1:] * kernel).mean()
    return dem - local_mean

def calculate_twi(slope, area):
    # Compute TWI based on the slope and contributing area
    twi = np.log(area / (np.tan(slope)))
    return twi

if dem_file and shapefile:
    with rasterio.open(dem_file) as src:
        dem = src.read(1)
        affine = src.transform

        # Compute slope using gradient
        x, y = np.gradient(dem, 1)
        slope = np.sqrt(x**2 + y**2)

        # Compute curvature (second derivative)
        curvature = np.gradient(slope, 1)

        # Compute TPI
        tpi = calculate_tpi(dem, window_size=3)

        # Compute catchment area approximation (as a proxy for flow accumulation)
        catchment_area = np.cumsum(dem)

        # Compute TWI
        twi = calculate_twi(slope, catchment_area)

        # Show the computed layers (Slope, Curvature, TPI, TWI)
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes[0, 0].set_title('Slope')
        show(slope, ax=axes[0, 0], transform=src.transform)

        axes[0, 1].set_title('Curvature')
        show(curvature, ax=axes[0, 1], transform=src.transform)

        axes[1, 0].set_title('TPI')
        show(tpi, ax=axes[1, 0], transform=src.transform)

        axes[1, 1].set_title('TWI')
        show(twi, ax=axes[1, 1], transform=src.transform)

        st.pyplot(fig)

        # Load shapefile and extract DEM values at points
        gdf = gpd.read_file(shapefile)
        st.write("Uploaded Point Shapefile:", gdf.head())

        # Extract values at points
        gdf['dem_value'] = [dem[src.index(x, y)] for x, y in zip(gdf.geometry.x, gdf.geometry.y)]
        gdf['slope_value'] = [slope[src.index(x, y)] for x, y in zip(gdf.geometry.x, gdf.geometry.y)]
        gdf['tpi_value'] = [tpi[src.index(x, y)] for x, y in zip(gdf.geometry.x, gdf.geometry.y)]
        gdf['twi_value'] = [twi[src.index(x, y)] for x, y in zip(gdf.geometry.x, gdf.geometry.y)]
        
        st.write("Extracted values at points:", gdf[['geometry', 'dem_value', 'slope_value', 'tpi_value', 'twi_value']])

        # Save the results to a new shapefile
        output_shapefile = "output_points_with_topographic_values.shp"
        gdf.to_file(output_shapefile)
        st.write(f"Saved extracted points with topographic values to: {output_shapefile}")


else:
    st.write("Please upload both a DEM and a point shapefile.")
