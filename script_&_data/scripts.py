import geopandas as gpd
from shapely.geometry import LineString, MultiLineString
from shapely.geometry import Point, MultiPoint
import os
import math
import warnings

# Ignore all warnings
warnings.filterwarnings("ignore")




def calculate_collocation_length(row):
    """get the collocation length of each pipeline"""
    return row['geometry'].length # calculate the leght of a line


def calculate_collocation_angle(pipeline_geometry, hight_tension_line):
    """
    Calculate the collocation angle between a pipeline and a high-tension power line.

    Args:
        pipeline_geometry (LineString): Shapely LineString representing the pipeline.
        hight_tension_line (LineString or MultiLineString): Shapely LineString or MultiLineString representing the high-tension line.

    Returns:
        float: The average collocation angle in degrees.
    """
    intersection_line = pipeline_geometry.intersection(hight_tension_line)
    
    if intersection_line.is_empty:
        return 0 # No intersection, return 0 angle
    
    elif isinstance(intersection_line, (Point, MultiPoint)):
        if isinstance(intersection_line, Point):
            intersection_coords = list(intersection_line.coords)[0]
            dx = pipeline_geometry.interpolate(pipeline_geometry.project(Point(intersection_coords))).x - intersection_coords[0]
            dy = pipeline_geometry.interpolate(pipeline_geometry.project(Point(intersection_coords))).y - intersection_coords[1]
            angle = math.degrees(math.atan2(dy, dx))
            return abs(angle)
        elif isinstance(intersection_line, MultiPoint):
            angles = []
            for point in intersection_line.geoms:
                intersection_coords = list(point.coords)[0]
                dx = pipeline_geometry.interpolate(pipeline_geometry.project(Point(intersection_coords))).x - intersection_coords[0]
                dy = pipeline_geometry.interpolate(pipeline_geometry.project(Point(intersection_coords))).y - intersection_coords[1]
                angle = math.degrees(math.atan2(dy, dx))
                angles.append(abs(angle))
            if angles:
                average_angle = sum(angles) / len(angles)
                return average_angle
            else:
                return 0  # No valid points, return 0 angle
    elif isinstance(intersection_line, (LineString, MultiLineString)):
        coords = list(intersection_line.coords)
        if len(coords) >= 2:
            dx = coords[-1][0] - coords[0][0]
            dy = coords[-1][1] - coords[0][1]
            angle = math.degrees(math.atan2(dy, dx))
            return abs(angle)
        else:
            return 0  # No valid line segments, return 0 angle
    else:
        return 0  # Other cases, return 0 angle

def calculate_separation_distance(pipeline_geometry: LineString, high_tension_geometry: LineString) -> float:
    """
    Calculate the separation distance between a pipeline and a high-tension power line.

    Args:
        pipeline_geometry (LineString): Shapely LineString representing the pipeline geometry.
        high_tension_geometry (LineString): Shapely LineString representing the high-tension power line geometry.

    Returns:
        float: The separation distance in feet.
    """
    # Calculate the minimum distance between the pipeline and high-tension line
    nearest_distance = pipeline_geometry.distance(high_tension_geometry)
    separation_distance= nearest_distance
    return separation_distance


def calculate_collocation_info(pipeline_gdf, high_tension_gdf, buffer_distance,target_crs):
    # Convert input data to the target CRS
    pipeline_gdf = pipeline_gdf.to_crs(target_crs)
    high_tension_gdf = high_tension_gdf.to_crs(target_crs)
    
    # Combine the high-tension line geometries into a single LineString
    merged_ht_geometry = high_tension_gdf.unary_union
    
    # Buffer around the merged high-tension power line
    buffer_polygon = merged_ht_geometry.buffer(buffer_distance)
    
    # Select pipeline features within the buffer
    pipeline_buffer_gdf = pipeline_gdf[pipeline_gdf.intersects(buffer_polygon)]
    
    # Calculate collocation length
    pipeline_buffer_gdf['col_length'] = pipeline_buffer_gdf.apply(calculate_collocation_length, axis=1)
    
    # Calculate collocation angle
    pipeline_buffer_gdf['col_angle'] = pipeline_buffer_gdf.apply(
        lambda row: calculate_collocation_angle(row['geometry'],merged_ht_geometry), axis=1)

    # Calculate separation distance
    pipeline_buffer_gdf['sep_dist'] = pipeline_buffer_gdf.apply(
        lambda row: calculate_separation_distance(row['geometry'],merged_ht_geometry), axis=1)
    
    # Select only the desired columns
    final_columns = ['PIPELINE_I', 'col_length', 'col_angle', 'sep_dist',
                     'geometry']
    pipeline_buffer_gdf = pipeline_buffer_gdf[final_columns]

    return pipeline_buffer_gdf


if __name__ == "__main__":
    # Load pipeline and high tension data
    path = 'données'
    pipeline_data = gpd.read_file(os.path.join(path,'pipeline1.shp'))
    high_tension_data = gpd.read_file(os.path.join(path,'H_T.shp'))
    
    # Set buffer distance, CRS, and thresholds
    buffer_distance = 800
    target_crs = "EPSG:3857"  # Define your target CRS here


    result = calculate_collocation_info(pipeline_data, high_tension_data, buffer_distance,target_crs)

    print(f"Number of pipeline features within {buffer_distance}m buffer of high tension line:", len(result))

    # Save result to a shapefile with GeoDataFrame export engine specified
    result.to_file('pipeline_collocation_result.shp', driver='ESRI Shapefile', encoding='utf-8')