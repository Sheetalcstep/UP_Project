r"""
Connect substations through the power-tower network from a CSV (batch).

Tower–tower graph edges:
- **0.5 km (500 m)**: any pair on land (default span).
- **1.2 km**: allowed only for that pair if the straight segment between the two
  towers **intersects** a water feature from ``water_bodies.shp`` (wider span only
  to bridge water).

Run from UP project folder: ``python code/Transmission_line_UP.py``

By default, reads next to this script: ``Demo.csv``,
``power_substation_UP_other_2.shp``, ``water_bodies.shp``. Override with
``UP_CSV``, ``UP_TOWER_SHP``, ``UP_OUT_SHP``, ``UP_FAILED_XLSX``, ``UP_WATER_SHP``,
and optionally ``UP_WATER_TOWER_MASK_BUFFER_M`` (meters, e.g. 700) to use a
smaller “who is near water” buffer for a faster bridge pass.
"""

from __future__ import annotations

import os
import time
import geopandas as gpd
import pandas as pd
import networkx as nx
from shapely.geometry import LineString
from shapely.ops import transform
from math import radians, cos, sin, asin, sqrt, atan2
from scipy.spatial import cKDTree
import numpy as np
import pyproj

# --- Tunable: tower–tower graph (per pair) -------------------------------------------
MAX_TOWER_LAND_KM = 0.5  # max edge when segment does not qualify as water span
MAX_TOWER_WATER_KM = 1.2  # max edge only if segment intersects a water body
MAX_SUBSTATION_TO_TOWER_KM = 1.0  # substation to nearest tower(s) (km)


def _script_dir() -> str:
    return os.path.dirname(os.path.abspath(__file__))


def _default_path(name: str) -> str:
    return os.path.join(_script_dir(), name)


# Defaults next to this file (override with env)
CSV_FILE = os.environ.get("UP_CSV", _default_path("excel/UP_power_map_line_length_csv_SC_DC.csv"))
TOWER_SHP = os.environ.get("UP_TOWER_SHP", _default_path("power_substation_UP_other_2.shp"))
WATER_SHP = os.environ.get("UP_WATER_SHP", _default_path("water_bodies.shp"))
OUTPUT_SHP = os.environ.get("UP_OUT_SHP", _default_path("Shapefile/UP_transmission_lines.shp"))
FAILED_XLSX = os.environ.get("UP_FAILED_XLSX", _default_path("excel/failed_central_UP.xlsx"))
WGS84 = "EPSG:4326"
# Meters. Candidates for the “j” tower in the bridge pass; smaller → faster (fewer j to
# loop). Default = full 1.2 km; set e.g. UP_WATER_TOWER_MASK_BUFFER_M=700 to shrink load.
WATER_TOWER_MASK_BUFFER_M = None
if os.environ.get("UP_WATER_TOWER_MASK_BUFFER_M", "").strip():
    try:
        WATER_TOWER_MASK_BUFFER_M = float(
            os.environ.get("UP_WATER_TOWER_MASK_BUFFER_M", "0")
        )
    except ValueError:
        WATER_TOWER_MASK_BUFFER_M = None


def _haversine_km_to_points(lj: float, tj: float, li: np.ndarray, la: np.ndarray) -> np.ndarray:
    """
    Great-circle distance (km) from one point to many. Vectorized, float64.
    """
    r = np.radians
    lj, tj = r(lj), r(tj)
    li = r(li.astype(float))
    la = r(la.astype(float))
    dlon = li - lj
    dlat = la - tj
    a = (
        np.sin(dlat * 0.5) ** 2
        + np.cos(tj) * np.cos(la) * (np.sin(dlon * 0.5) ** 2)
    )
    a = np.clip(a, 0.0, 1.0)
    c = 2.0 * np.arcsin(np.sqrt(a))
    return 6371.0 * c


def _line_bbox_overlaps_extent(
    lon1: float, lat1: float, lon2: float, lat2: float, extent: tuple[float, float, float, float]
) -> bool:
    """
    If false, a segment between the two WGS84 points cannot intersect *any* geometry
    with bbox ``extent = (minx, miny, maxx, maxy)`` (e.g. water layer).
    """
    wminx, wminy, wmaxx, wmaxy = extent
    Lminx, Lminy = min(lon1, lon2), min(lat1, lat2)
    Lmaxx, Lmaxy = max(lon1, lon2), max(lat1, lat2)
    if Lmaxx < wminx or Lminx > wmaxx or Lmaxy < wminy or Lminy > wmaxy:
        return False
    return True


def haversine_distance(lon1, lat1, lon2, lat2):
    """
    Great-circle distance in km between two WGS84 points in decimal degrees.
    """
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    r = 6371.0
    return c * r


def _degrees(rad):
    return rad * 180.0 / 3.141592653589793


def create_simple_offset(line, offset_meters=10):
    """
    Simple perpendicular offset for small distances (fallback).
    """
    coords = list(line.coords)
    if len(coords) < 2:
        return line
    lon1, lat1 = coords[0]
    lon2, lat2 = coords[-1]
    lat1_rad, lon1_rad = radians(lat1), radians(lon1)
    lat2_rad, lon2_rad = radians(lat2), radians(lon2)
    dlon = lon2_rad - lon1_rad
    bearing = atan2(
        sin(dlon) * cos(lat2_rad),
        cos(lat1_rad) * sin(lat2_rad) - sin(lat1_rad) * cos(lat2_rad) * cos(dlon),
    )
    perp_bearing = bearing + radians(90)
    offset_deg = offset_meters / 111000.0
    offset_coords = []
    for lon, lat in coords:
        lat_rad = radians(lat)
        dlat_offset = offset_deg * cos(perp_bearing)
        dlon_offset = offset_deg * sin(perp_bearing) / cos(lat_rad)
        offset_coords.append((lon + _degrees(dlon_offset), lat + _degrees(dlat_offset)))
    return LineString(offset_coords)


def create_parallel_line(line, offset_meters=10):
    """
    Offset a LineString by ``offset_meters`` (UTM, then back to WGS84).
    """
    centroid = line.centroid
    lon, lat = centroid.x, centroid.y
    utm_zone = int((lon + 180) / 6) + 1
    if lat >= 0:
        utm_crs = f"EPSG:326{utm_zone:02d}"
        
    else:
        utm_crs = f"EPSG:327{utm_zone:02d}"
    try:
        project_to_utm = pyproj.Transformer.from_crs(
            "EPSG:4326", utm_crs, always_xy=True
        ).transform
        project_to_wgs = pyproj.Transformer.from_crs(
            utm_crs, "EPSG:4326", always_xy=True
        ).transform
        line_utm = transform(project_to_utm, line)
        offset_line_utm = line_utm.parallel_offset(offset_meters, "right")
        return transform(project_to_wgs, offset_line_utm)
    except Exception as e:
        print(f"Warning: UTM projection failed, using simple offset: {e}")
        return create_simple_offset(line, offset_meters)


def load_water_bodies(path: str) -> gpd.GeoDataFrame:
    """
    Load water features and reproject to WGS84 (lon, lat) for use with tower segments.
    """
    gdf = gpd.read_file(path)
    gdf = gdf[~gdf.geometry.isna()]
    gdf = gdf[~gdf.geometry.is_empty]
    if gdf.crs is None or not gdf.crs:
        gdf = gdf.set_crs(WGS84)
    else:
        gdf = gdf.to_crs(WGS84)
    gdf = gdf.reset_index(drop=True)
    return gdf


def _build_water_strtree(
    water: gpd.GeoDataFrame | None,
) -> tuple[object, list, gpd.GeoDataFrame] | None:
    """
    Build a Shapely STRtree (plus geometry list) once for many segment checks.
    Falls back to (None, [], water) if STRtree is unavailable; queries use sjoin.
    """
    if water is None or water.empty:
        return None
    geoms = [
        g
        for g in water.geometry
        if g is not None
        and not (hasattr(g, "is_empty") and g.is_empty)
    ]
    if not geoms:
        return None
    try:
        from shapely.strtree import STRtree
    except ImportError:  # pragma: no cover
        STRtree = None
    if STRtree is not None:
        try:
            return (STRtree(geoms), geoms, water)
        except Exception:  # pragma: no cover
            pass
    return (None, geoms, water)


def _segment_crosses_water(
    lon1: float, lat1: float, lon2: float, lat2: float,
    w_index: tuple[object, list, gpd.GeoDataFrame] | None,
) -> bool:
    """True if the tower–tower great-circle hop (WGS84 lon, lat) intersects water features."""
    if w_index is None:
        return False
    tree, geoms, water_gdf = w_index[0], w_index[1], w_index[2]
    line = LineString([(lon1, lat1), (lon2, lat2)])
    if not line.is_valid or line.is_empty or line.length == 0:
        return False
    if tree is not None:
        try:
            idxs = tree.query(line, predicate="intersects")
        except TypeError:
            try:
                idxs = tree.query(line)
            except Exception:
                idxs = np.array([], dtype=int)
        for k in np.atleast_1d(idxs).ravel().astype(int):
            g = geoms[int(k)]
            if g is not None and not (getattr(g, "is_empty", True)) and g.intersects(line):
                return True
        return False
    g_line = gpd.GeoDataFrame(geometry=[line], crs=water_gdf.crs)
    j = gpd.sjoin(
        g_line, water_gdf, how="inner", predicate="intersects"
    )
    return not j.empty


def _utm_crs_for_mean_lonlat(lon: float, lat: float) -> str:
    """Northern UTM zone from a sample lon/lat (WGS84)."""
    z = int((float(lon) + 180.0) / 6.0) + 1
    if float(lat) >= 0:
        return f"EPSG:326{min(60, max(1, z)):02d}"
    return f"EPSG:327{min(60, max(1, z)):02d}"


def _tower_pos_near_water_mask(
    lons: np.ndarray,
    lats: np.ndarray,
    water: gpd.GeoDataFrame,
    buffer_m: float,
) -> np.ndarray:
    """
    True for tower row positions (0..n-1) within ``buffer_m`` of any water geometry
    (in projected meters). Used to limit expensive 0.5–1.2 km pair checks to plausible
    river-bank candidates. Uses a spatial join (no union of all water), so it scales
    to large water layers.
    """
    n = len(lons)
    if n == 0 or water is None or water.empty:
        return np.zeros(0, dtype=bool)
    m_lon, m_lat = float(np.mean(lons)), float(np.mean(lats))
    utm = _utm_crs_for_mean_lonlat(m_lon, m_lat)
    try:
        gxy = gpd.points_from_xy(lons, lats, crs=WGS84)
    except Exception:  # pragma: no cover
        from shapely.geometry import Point
        gxy = [Point(float(x), float(y)) for x, y in zip(lons, lats)]
    pts = gpd.GeoDataFrame(
        {"_pid": np.arange(n, dtype=np.int64)},
        geometry=gxy,
        crs=WGS84,
    ).to_crs(utm)
    w_utm = water.to_crs(utm)
    w_buf = w_utm.copy()
    w_buf["geometry"] = w_buf["geometry"].buffer(float(buffer_m))
    joined = gpd.sjoin(
        pts, w_buf[["geometry"]], how="inner", predicate="intersects"
    )
    out = np.zeros(n, dtype=bool)
    if len(joined) == 0:
        return out
    for pid in np.unique(joined["_pid"].to_numpy()):
        p = int(pid)
        if 0 <= p < n:
            out[p] = True
    return out


def build_tower_network(
    towers_gdf,
    water: gpd.GeoDataFrame | None = None,
    max_tower_land_km: float = MAX_TOWER_LAND_KM,
    max_tower_water_km: float = MAX_TOWER_WATER_KM,
):
    """
    Build an undirected graph. Edges up to ``max_tower_land_km`` (km) in general; between
    ``max_tower_land_km`` and ``max_tower_water_km`` only if the hop intersects a water
    feature from the loaded water layer.

    For speed, land edges use a small ``query_pairs`` radius; water-only spans are
    checked only for towers within a ~``max_tower_water_km`` buffer of water, then
    confirmed with a line–water intersection (not a global 1.2 km pair search).
    """
    w_index = _build_water_strtree(water)
    g = nx.Graph()
    tower_coords: dict[str, tuple[float, float]] = {}
    lons = towers_gdf["Long"].values
    lats = towers_gdf["Lat"].values
    n = len(lons)
    for idx, row in towers_gdf.iterrows():
        node_id = f"tower_{idx}"
        lon = float(row["Long"])
        lat = float(row["Lat"])
        tower_coords[node_id] = (lon, lat)
        g.add_node(node_id, lon=lon, lat=lat, idx=idx)
    avg_lat = float(np.mean(lats))
    lat_deg_per_km = 1.0 / 111.0
    lon_deg_per_km = 1.0 / (111.0 * np.cos(np.radians(avg_lat)))
    scale = max(lat_deg_per_km, lon_deg_per_km)
    coords_array = np.column_stack([lons, lats])
    tree = cKDTree(coords_array)
    # --- Phase 1: all pairs within land span (small ring → few pairs) ----------------
    r_land = max_tower_land_km * scale * 1.5
    pairs_land = tree.query_pairs(r_land, output_type="ndarray")
    land_edges = 0
    for i, j in pairs_land:
        idx1 = towers_gdf.index[i]
        idx2 = towers_gdf.index[j]
        node1 = f"tower_{idx1}"
        node2 = f"tower_{idx2}"
        a1, b1, a2, b2 = lons[i], lats[i], lons[j], lats[j]
        d = haversine_distance(a1, b1, a2, b2)
        if d <= max_tower_land_km:
            g.add_edge(node1, node2, weight=d)
            land_edges += 1
    # --- Phase 2: 0.5–1.2 km only for towers near water; verify segment vs water ----
    water_edges = 0
    if w_index is not None and water is not None and not water.empty:
        mask_buf = (
            WATER_TOWER_MASK_BUFFER_M
            if WATER_TOWER_MASK_BUFFER_M is not None
            else max_tower_water_km * 1000.0
        )
        print(
            f"  Marking towers near water (buffer {mask_buf/1000.0} km)…",
            flush=True,
        )
        near = _tower_pos_near_water_mask(
            lons, lats, water, buffer_m=mask_buf
        )
        njp = int(np.count_nonzero(near))
        print(
            f"  {njp} / {n} tower positions in water buffer; checking bridge hops (i<j, vector dist)…",
            flush=True,
        )
        r_w = max_tower_water_km * scale * 1.5
        w_ext = water.total_bounds
        idx_to_label = towers_gdf.index.to_numpy()
        near_list = np.flatnonzero(near)
        t_bridge = time.time()
        for step, j in enumerate(near_list, start=1):
            if (step > 0 and (step % 10000) == 0) or step == njp or step == 1:
                elapsed = time.time() - t_bridge
                rate = step / max(elapsed, 0.1)
                left = (njp - step) / max(rate, 1e-6)
                print(
                    f"  … bridge {step}/{njp} (~{left/60.0:.1f} min left, "
                    f"elapsed {int(elapsed)}s)",
                    flush=True,
                )
            cands = list(tree.query_ball_point(coords_array[j], r_w))
            if not cands:
                continue
            ix = np.fromiter((i for i in cands if i < j), dtype=np.int64)
            if ix.size == 0:
                continue
            lj, tj = lons[j], lats[j]
            d = _haversine_km_to_points(lj, tj, lons[ix], lats[ix])
            m = (d > max_tower_land_km) & (d <= max_tower_water_km)
            for k, w in zip(ix[m], d[m]):
                lo1, la1, lo2, la2 = lons[k], lats[k], lons[j], lats[j]
                if not _line_bbox_overlaps_extent(lo1, la1, lo2, la2, w_ext):
                    continue
                if not _segment_crosses_water(lo1, la1, lo2, la2, w_index):
                    continue
                p, q = int(k), int(j)
                lab1, lab2 = idx_to_label[p], idx_to_label[q]
                n1, n2 = f"tower_{lab1}", f"tower_{lab2}"
                g.add_edge(n1, n2, weight=float(w))
                water_edges += 1
    print(
        f"  Graph: {land_edges} land edges (≤{max_tower_land_km} km), "
        f"{water_edges} water bridge edges",
        flush=True,
    )
    return g, tower_coords


def find_nearest_towers(
    substation_lon,
    substation_lat,
    towers_gdf,
    max_distance: float = 2.0,
):
    lons = towers_gdf["Long"].values
    lats = towers_gdf["Lat"].values
    avg_lat = float(np.mean(lats))
    lat_deg_per_km = 1.0 / 111.0
    lon_deg_per_km = 1.0 / (111.0 * np.cos(np.radians(avg_lat)))
    max_deg = max_distance * max(lat_deg_per_km, lon_deg_per_km)
    tree = cKDTree(np.column_stack([lons, lats]))
    query_point = np.array([substation_lon, substation_lat])
    out = []
    for idx_pos in tree.query_ball_point(query_point, max_deg * 1.5):
        idx = towers_gdf.index[idx_pos]
        tlon, tlat = lons[idx_pos], lats[idx_pos]
        d = haversine_distance(substation_lon, substation_lat, tlon, tlat)
        if d <= max_distance:
            out.append(f"tower_{idx}")
    return out


def connect_substations_from_coords(
    from_lon,
    from_lat,
    to_lon,
    to_lat,
    towers_gdf,
    g,
    tower_coords,
    max_substation_distance: float = 2.0,
):
    from_towers = find_nearest_towers(
        from_lon, from_lat, towers_gdf, max_substation_distance
    )
    to_towers = find_nearest_towers(
        to_lon, to_lat, towers_gdf, max_substation_distance
    )
    if not from_towers:
        return None, None, None, "No towers near from substation"
    if not to_towers:
        return None, None, None, "No towers near to substation"
    g_conn = g.copy()
    g_conn.add_node("from_sub", lon=from_lon, lat=from_lat)
    g_conn.add_node("to_sub", lon=to_lon, lat=to_lat)
    for tower_node in from_towers:
        tlon, tlat = tower_coords[tower_node]
        g_conn.add_edge(
            "from_sub", tower_node, weight=haversine_distance(
                from_lon, from_lat, tlon, tlat
            )
        )
    for tower_node in to_towers:
        tlon, tlat = tower_coords[tower_node]
        g_conn.add_edge(
            "to_sub", tower_node, weight=haversine_distance(
                to_lon, to_lat, tlon, tlat
            )
        )
    try:
        path = nx.shortest_path(
            g_conn, "from_sub", "to_sub", weight="weight"
        )
        plen = nx.shortest_path_length(
            g_conn, "from_sub", "to_sub", weight="weight"
        )
        path_coords = []
        for node in path:
            if node == "from_sub":
                path_coords.append((from_lon, from_lat))
            elif node == "to_sub":
                path_coords.append((to_lon, to_lat))
            else:
                path_coords.append(tower_coords[node])
        path_line = LineString(path_coords)
        return path_coords, path_line, plen, None
    except nx.NetworkXNoPath:
        return None, None, None, "No path found through tower network"


def main():
    max_substation_distance = MAX_SUBSTATION_TO_TOWER_KM
    csv_file = CSV_FILE
    tower_shp = TOWER_SHP
    water_shp = WATER_SHP
    output_shp = OUTPUT_SHP
    failed_connections_file = FAILED_XLSX

    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"CSV file not found: {csv_file}")
    if not os.path.exists(tower_shp):
        raise FileNotFoundError(f"Tower shapefile not found: {tower_shp}")
    if not os.path.exists(water_shp):
        raise FileNotFoundError(
            f"Water shapefile not found: {water_shp}. Set UP_WATER_SHP or add water_bodies.shp."
        )

    print("=" * 80)
    print("BATCH PROCESSING (500 m land / 1.2 km water-hops per pair)")
    print("=" * 80)

    print(f"\nLoading CSV file: {csv_file}")
    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} connections from CSV")

    print(f"\nLoading tower data: {tower_shp}")
    towers_gdf = gpd.read_file(tower_shp)
    print(f"Loaded {len(towers_gdf)} towers")

    if "Lat" not in towers_gdf.columns or "Long" not in towers_gdf.columns:
        print("Tower columns:", list(towers_gdf.columns))
        raise ValueError("Tower shapefile must have 'Lat' and 'Long' columns")

    required_cols = [
        "Fr_SS_Lat", "Fr_SS_Long", "To_SS_Lat", "To_SS_Long", "No_of_ckt"
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print("CSV columns:", list(df.columns))
        raise ValueError(f"CSV file must have these columns: {missing}")

    print(f"\nLoading water: {water_shp}")
    water_gdf = load_water_bodies(water_shp)
    print(
        f"  {len(water_gdf)} water features in WGS84; "
        f"edges: ≤{MAX_TOWER_LAND_KM} km land, "
        f"≤{MAX_TOWER_WATER_KM} km if segment hits water"
    )
    print("\nBuilding tower network...")
    g, tower_coords = build_tower_network(
        towers_gdf,
        water=water_gdf,
        max_tower_land_km=MAX_TOWER_LAND_KM,
        max_tower_water_km=MAX_TOWER_WATER_KM,
    )
    print(
        f"Network built: {g.number_of_nodes()} nodes, {g.number_of_edges()} edges"
    )

    print(f"\nProcessing {len(df)} connections...")
    print("=" * 80)

    all_lines: list[dict] = []
    failed_connections: list[dict] = []

    for idx, row in df.iterrows():
        try:
            from_lat = row["Fr_SS_Lat"]
            from_lon = row["Fr_SS_Long"]
            to_lat = row["To_SS_Lat"]
            to_lon = row["To_SS_Long"]
            circuits = str(row["No_of_ckt"]).strip().upper()

            if (
                pd.isna(from_lat)
                or pd.isna(from_lon)
                or pd.isna(to_lat)
                or pd.isna(to_lon)
            ):
                error_msg = f"Row {idx + 1}: Missing coordinates"
                print(f"! {error_msg}")
                failed_connections.append({**row.to_dict(), "Error": error_msg})
                continue

            if str(from_lon) == "#N/A" or str(to_lon) == "#N/A":
                error_msg = f"Row {idx + 1}: Invalid coordinate (#N/A)"
                print(f"! {error_msg}")
                failed_connections.append({**row.to_dict(), "Error": error_msg})
                continue

            from_lat, from_lon = float(from_lat), float(from_lon)
            to_lat, to_lon = float(to_lat), float(to_lon)

            print(
                f"\nProcessing row {idx + 1}/{len(df)}: "
                f"{row.get('Fr_SS', 'N/A')} -> {row.get('To_SS', 'N/A')}"
            )

            path_coords, path_line, total_distance, err = connect_substations_from_coords(
                from_lon,
                from_lat,
                to_lon,
                to_lat,
                towers_gdf,
                g,
                tower_coords,
                max_substation_distance,
            )

            if err:
                print(f"  Failed: {err}")
                failed_connections.append({**row.to_dict(), "Error": err})
                continue

            print(
                f"  Path: {total_distance:.2f} km, {len(path_coords)} points"
            )

            if circuits == "S/C":
                line_data = row.to_dict()
                line_data["distance_km"] = total_distance
                line_data["num_towers"] = len(path_coords) - 2
                all_lines.append({"geometry": path_line, **line_data})
            elif circuits == "D/C":
                line1_data = row.to_dict()
                line1_data["distance_km"] = total_distance
                line1_data["num_towers"] = len(path_coords) - 2
                line1_data["circuit"] = "Line1"
                all_lines.append({"geometry": path_line, **line1_data})
                try:
                    offset_line = create_parallel_line(path_line, offset_meters=10)
                    line2_data = row.to_dict()
                    line2_data["distance_km"] = total_distance
                    line2_data["num_towers"] = len(path_coords) - 2
                    line2_data["circuit"] = "Line2"
                    all_lines.append(
                        {"geometry": offset_line, **line2_data}
                    )
                    print("  Parallel line (10 m offset)")
                except Exception as e:
                    print(
                        f"  Warning: Could not create parallel line: {e}"
                    )
                    line2_data = row.to_dict()
                    line2_data["distance_km"] = total_distance
                    line2_data["num_towers"] = len(path_coords) - 2
                    line2_data["circuit"] = "Line2"
                    all_lines.append({"geometry": path_line, **line2_data})
            else:
                print(
                    f"  Warning: Unknown circuit type {circuits!r}, using single line"
                )
                line_data = row.to_dict()
                line_data["distance_km"] = total_distance
                line_data["num_towers"] = len(path_coords) - 2
                all_lines.append({"geometry": path_line, **line_data})
        except Exception as e:
            error_msg = f"Row {idx + 1}: {e!s}"
            print(f"  Error: {error_msg}")
            failed_connections.append(
                {**row.to_dict(), "Error": error_msg}
            )
            continue

    print(f"\n{'='*80}")
    print("Creating output shapefile...")
    print(f"Successfully processed: {len(all_lines)} lines")
    print(f"Failed connections: {len(failed_connections)}")

    if all_lines:
        output_gdf = gpd.GeoDataFrame(all_lines, crs=towers_gdf.crs)
        base = output_shp.replace(".shp", "")
        for ext in (".shp", ".shx", ".dbf", ".prj", ".cpg", ".sbn", ".sbx"):
            p = base + ext
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        output_gdf.to_file(output_shp)
        print(f"\nLines saved: {output_shp} ({len(output_gdf)} features)")
    else:
        print("\nNo lines were created.")

    if failed_connections:
        failed_df = pd.DataFrame(failed_connections)
        failed_df.to_excel(failed_connections_file, index=False)
        print(
            f"\nFailed connections: {failed_connections_file} "
            f"({len(failed_connections)} rows)"
        )
    else:
        print("\nAll connections OK.")

    print(f"\n{'='*80}\nDone.\n{'='*80}")


if __name__ == "__main__":
    main()
