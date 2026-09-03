r"""
Route failed substation connections through the power-tower network -- with a
flagged "bridge" fallback for genuine gaps in the tower data.

This is a new script, separate from route_failed_connections_v2.py (left
untouched) and Transmission_line_UP.py (also untouched). Everything about the
normal routing is identical to v2:

Tower-tower graph edges:
- 1 km flat. No water-crossing exception.

Substation-tower edges:
- Expanding rings: 0-1 km, 1-2 km, 2-3 km, ... up to MAX_SUBSTATION_SEARCH_KM.
  Connects to every tower found in the first non-empty ring.

What's new -- the bridge fallback:
Some rows fail not because a substation has no nearby tower, but because the
towers it DOES find sit in their own small, disconnected island -- cut off
from the rest of the tower network by a real gap in the source data (missing
towers in that stretch). Widening the 1 km tower-to-tower rule for everyone
to fix this is risky (that's what caused the zigzag problems before), so
instead:

1. Try the normal path first, exactly like v2.
2. If that fails, find the two tower "islands" (network components) the from
   and to substations actually landed in, and measure the real gap between
   them -- the closest distance between any tower in one island and any
   tower in the other.
3. If that gap is <= BRIDGE_MAX_KM, add ONE straight bridge segment across
   just that gap, connect the islands, and re-route. The output row is
   tagged with a `bridge_km` value so you can see exactly which lines
   contain an auto-bridged stretch and go check it on the map.
4. If the gap is bigger than BRIDGE_MAX_KM, the row is left as failed for
   manual review, same as before -- but the error message now also reports
   how big the real gap actually is, so you know whether it's close or way
   off.

BRIDGE_MAX_KM defaults to 8 km (comfortably above the 5.3 km gap found for
the SHAHJAHAN PUR(PG) -> MALLAWAN case). Change it near the top of this file.

Run from the UP project folder:
    python route_failed_connections_v3.py

Reads (override with UP3_FAILED_XLSX / UP3_TOWER_SHP / UP3_OUT_SHP /
UP3_REMAINING_XLSX / UP3_ROW_LIMIT env vars):
    excel/failed_UP_3_remaining.xlsx   (currently the same list v2 is working from)
    power_substation_UP_other_2.shp

Writes:
    Shapefile/failed_lines_UP_bridged.shp
    excel/failed_UP_bridged_remaining.xlsx
    logs/route_v3_<timestamp>.log

ROW_LIMIT near the top controls how many rows get processed -- set to 10 for
a quick test batch. Set it to None (or run with UP3_ROW_LIMIT=0) for all rows.
"""

from __future__ import annotations

import datetime
import os
import sys
import time
from pathlib import Path

import geopandas as gpd
import pandas as pd
import networkx as nx
from shapely.geometry import LineString
from shapely.ops import transform
from math import radians, cos, sin, asin, sqrt, atan2
from scipy.spatial import cKDTree
import numpy as np
import pyproj

# --- Tunable distances ------------------------------------------------------------------
MAX_TOWER_KM = 1.0               # tower-to-tower graph edge, flat (no water exception)
SUBSTATION_RING_STEP_KM = 1.0    # substation-to-tower search ring width
MAX_SUBSTATION_SEARCH_KM = 5.0  # substation-to-tower search cap (0-1, 1-2, ... up to this)
BRIDGE_MAX_KM = 6.0              # max gap the fallback is allowed to auto-bridge

# Test-batch size: only process the first N rows of the excel file so you can check
# the output before committing to a full run. Set to None (or override with the
# UP3_ROW_LIMIT env var, e.g. UP3_ROW_LIMIT=0 for "all rows") once you're happy.
ROW_LIMIT = 12


WGS84 = "EPSG:4326"


def _script_dir() -> Path:
    return Path(__file__).resolve().parent


def _default_path(name: str) -> str:
    return str(_script_dir() / name)


FAILED_XLSX = os.environ.get("UP3_FAILED_XLSX", _default_path("excel/failed/New folder/UP_41_lines_remaining.xlsx"))
TOWER_SHP = os.environ.get("UP3_TOWER_SHP", _default_path("power_substation_UP_other_2.shp"))
OUTPUT_SHP = os.environ.get("UP3_OUT_SHP", _default_path("Shapefile/UP_41_lines_remaining.shp"))
REMAINING_FAILED_XLSX = os.environ.get(
    "UP3_REMAINING_XLSX", _default_path("excel/failed/New folder/UP_17_lines_remaining_2.xlsx")
)

if os.environ.get("UP3_ROW_LIMIT", "").strip():
    _env_limit = int(os.environ["UP3_ROW_LIMIT"])
    ROW_LIMIT = None if _env_limit <= 0 else _env_limit


class Tee:
    """Writes to two streams at once (terminal + log file), flushing immediately."""

    def __init__(self, *streams):
        self._streams = streams

    def write(self, data):
        for s in self._streams:
            s.write(data)
            s.flush()

    def flush(self):
        for s in self._streams:
            s.flush()


# --------------------------------------------------------------------------------------
# Geometry helpers
# --------------------------------------------------------------------------------------

def haversine_distance(lon1, lat1, lon2, lat2):
    """Great-circle distance in km between two WGS84 points in decimal degrees."""
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * asin(sqrt(a))
    return 6371.0 * c


def _haversine_vec(lon1, lat1, lon2_arr, lat2_arr):
    """Great-circle distance (km) from one point to an array of points."""
    lon1r, lat1r = radians(lon1), radians(lat1)
    lon2r = np.radians(lon2_arr)
    lat2r = np.radians(lat2_arr)
    dlon = lon2r - lon1r
    dlat = lat2r - lat1r
    a = np.sin(dlat / 2) ** 2 + cos(lat1r) * np.cos(lat2r) * (np.sin(dlon / 2) ** 2)
    a = np.clip(a, 0.0, 1.0)
    return 6371.0 * 2.0 * np.arcsin(np.sqrt(a))


def _degrees(rad):
    return rad * 180.0 / 3.141592653589793


def create_simple_offset(line, offset_meters=10):
    """Perpendicular offset fallback if UTM projection fails."""
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
    """Offset a LineString by offset_meters (project to UTM, offset, project back)."""
    centroid = line.centroid
    lon, lat = centroid.x, centroid.y
    utm_zone = int((lon + 180) / 6) + 1
    utm_crs = f"EPSG:326{utm_zone:02d}" if lat >= 0 else f"EPSG:327{utm_zone:02d}"
    try:
        project_to_utm = pyproj.Transformer.from_crs(WGS84, utm_crs, always_xy=True).transform
        project_to_wgs = pyproj.Transformer.from_crs(utm_crs, WGS84, always_xy=True).transform
        line_utm = transform(project_to_utm, line)
        offset_line_utm = line_utm.parallel_offset(offset_meters, "right")
        return transform(project_to_wgs, offset_line_utm)
    except Exception as e:
        print(f"  Warning: UTM projection failed, using simple offset: {e}")
        return create_simple_offset(line, offset_meters)


# --------------------------------------------------------------------------------------
# Tower network
# --------------------------------------------------------------------------------------

def build_tower_network(towers_gdf, max_tower_km=MAX_TOWER_KM):
    """
    Undirected graph. Every pair of towers within max_tower_km (flat, no
    water-crossing exception) gets an edge weighted by real distance.
    """
    g = nx.Graph()
    tower_coords: dict[str, tuple[float, float]] = {}
    lons = towers_gdf["Long"].values
    lats = towers_gdf["Lat"].values

    for idx, row in towers_gdf.iterrows():
        node_id = f"tower_{idx}"
        lon, lat = float(row["Long"]), float(row["Lat"])
        tower_coords[node_id] = (lon, lat)
        g.add_node(node_id, lon=lon, lat=lat, idx=idx)

    avg_lat = float(np.mean(lats))
    lat_deg_per_km = 1.0 / 111.0
    lon_deg_per_km = 1.0 / (111.0 * np.cos(np.radians(avg_lat)))
    scale = max(lat_deg_per_km, lon_deg_per_km)
    coords_array = np.column_stack([lons, lats])
    tree = cKDTree(coords_array)

    r = max_tower_km * scale * 1.5
    pairs = tree.query_pairs(r, output_type="ndarray")
    edges = 0
    for i, j in pairs:
        idx1, idx2 = towers_gdf.index[i], towers_gdf.index[j]
        node1, node2 = f"tower_{idx1}", f"tower_{idx2}"
        d = haversine_distance(lons[i], lats[i], lons[j], lats[j])
        if d <= max_tower_km:
            g.add_edge(node1, node2, weight=d)
            edges += 1

    print(f"  Graph: {edges} edges (<= {max_tower_km} km)", flush=True)
    return g, tower_coords


def find_towers_within(substation_lon, substation_lat, towers_gdf, tower_tree, tower_lons, tower_lats,
                        max_distance):
    """All towers within max_distance km of the substation."""
    avg_lat = float(np.mean(tower_lats))
    lat_deg_per_km = 1.0 / 111.0
    lon_deg_per_km = 1.0 / (111.0 * np.cos(np.radians(avg_lat)))
    max_deg = max_distance * max(lat_deg_per_km, lon_deg_per_km)
    query_point = np.array([substation_lon, substation_lat])
    out = []
    for idx_pos in tower_tree.query_ball_point(query_point, max_deg * 1.5):
        idx = towers_gdf.index[idx_pos]
        tlon, tlat = tower_lons[idx_pos], tower_lats[idx_pos]
        d = haversine_distance(substation_lon, substation_lat, tlon, tlat)
        if d <= max_distance:
            out.append(f"tower_{idx}")
    return out


def find_towers_incremental(substation_lon, substation_lat, towers_gdf, tower_tree, tower_lons, tower_lats,
                             ring_step_km=SUBSTATION_RING_STEP_KM, max_search_km=MAX_SUBSTATION_SEARCH_KM):
    """
    Expand the substation-to-tower search in rings: 0-1 km, 1-2 km, 2-3 km, ...
    up to max_search_km. As soon as a ring's outer edge contains at least one
    tower, return every tower within that radius (not just the closest one)
    plus the radius it took to find them. Returns ([], None) if nothing is
    found within max_search_km at all.
    """
    ring_end = ring_step_km
    while ring_end <= max_search_km + 1e-9:
        towers = find_towers_within(
            substation_lon, substation_lat, towers_gdf, tower_tree, tower_lons, tower_lats, ring_end
        )
        if towers:
            return towers, ring_end
        ring_end += ring_step_km
    return [], None


# --------------------------------------------------------------------------------------
# Bridge fallback: find the two disconnected islands and measure the real gap
# --------------------------------------------------------------------------------------

def unique_components_for(g, node_list):
    """
    Return the distinct connected components (as frozensets of node ids) that
    the given towers belong to, deduplicated (several towers often share the
    same component).
    """
    seen = set()
    components = []
    for node in node_list:
        if node in seen:
            continue
        comp = nx.node_connected_component(g, node)
        seen.update(comp)
        components.append(frozenset(comp))
    return components


def nearest_gap_between(tower_coords, comp_a, comp_b):
    """
    Closest real distance (km) between any tower in comp_a and any tower in
    comp_b, plus the two tower ids that mark that gap. Iterates over the
    smaller component (vectorized with numpy against the larger one) so this
    stays fast even when one side is a large, well-connected network and the
    other is a small isolated island -- the common real-world case here.
    """
    if len(comp_a) > len(comp_b):
        comp_a, comp_b = comp_b, comp_a
    b_list = list(comp_b)
    b_lons = np.array([tower_coords[n][0] for n in b_list])
    b_lats = np.array([tower_coords[n][1] for n in b_list])

    best_dist = float("inf")
    best_a, best_b = None, None
    for a_node in comp_a:
        a_lon, a_lat = tower_coords[a_node]
        d = _haversine_vec(a_lon, a_lat, b_lons, b_lats)
        j = int(np.argmin(d))
        if d[j] < best_dist:
            best_dist = float(d[j])
            best_a, best_b = a_node, b_list[j]
    return best_dist, best_a, best_b


def connect_substations_from_coords(from_lon, from_lat, to_lon, to_lat,
                                     towers_gdf, tower_tree, tower_lons, tower_lats,
                                     g, tower_coords,
                                     ring_step_km=SUBSTATION_RING_STEP_KM,
                                     max_search_km=MAX_SUBSTATION_SEARCH_KM,
                                     bridge_max_km=BRIDGE_MAX_KM):
    """
    Returns: path_coords, path_line, total_distance, error, from_ring, to_ring, bridge_km
    bridge_km is None unless the fallback bridge was actually used.
    """
    from_towers, from_ring = find_towers_incremental(
        from_lon, from_lat, towers_gdf, tower_tree, tower_lons, tower_lats, ring_step_km, max_search_km
    )
    to_towers, to_ring = find_towers_incremental(
        to_lon, to_lat, towers_gdf, tower_tree, tower_lons, tower_lats, ring_step_km, max_search_km
    )
    if not from_towers:
        return None, None, None, f"No tower found within {max_search_km}km of from substation", None, None, None
    if not to_towers:
        return None, None, None, f"No tower found within {max_search_km}km of to substation", None, None, None

    g_conn = g.copy()
    g_conn.add_node("from_sub", lon=from_lon, lat=from_lat)
    g_conn.add_node("to_sub", lon=to_lon, lat=to_lat)
    for tower_node in from_towers:
        tlon, tlat = tower_coords[tower_node]
        g_conn.add_edge("from_sub", tower_node, weight=haversine_distance(from_lon, from_lat, tlon, tlat))
    for tower_node in to_towers:
        tlon, tlat = tower_coords[tower_node]
        g_conn.add_edge("to_sub", tower_node, weight=haversine_distance(to_lon, to_lat, tlon, tlat))

    bridge_km = None

    try:
        path = nx.shortest_path(g_conn, "from_sub", "to_sub", weight="weight")
    except nx.NetworkXNoPath:
        # Normal path failed -- find the real gap between the two islands and
        # try to bridge it, instead of giving up immediately.
        from_components = unique_components_for(g, from_towers)
        to_components = unique_components_for(g, to_towers)

        best_gap = (float("inf"), None, None)
        for comp_a in from_components:
            for comp_b in to_components:
                if comp_a & comp_b:
                    continue  # already the same island -- shouldn't happen if shortest_path failed
                gap_km, node_a, node_b = nearest_gap_between(tower_coords, comp_a, comp_b)
                if gap_km < best_gap[0]:
                    best_gap = (gap_km, node_a, node_b)

        gap_km, node_a, node_b = best_gap
        if node_a is None or gap_km > bridge_max_km:
            gap_desc = f"{gap_km:.2f} km" if node_a is not None else "unknown"
            return (None, None, None,
                    f"No path found through tower network (nearest gap between clusters: {gap_desc}, "
                    f"exceeds bridge cap of {bridge_max_km}km)",
                    from_ring, to_ring, None)

        # Bridge the gap with one extra edge and retry.
        g_conn.add_edge(node_a, node_b, weight=gap_km)
        bridge_km = gap_km
        try:
            path = nx.shortest_path(g_conn, "from_sub", "to_sub", weight="weight")
        except nx.NetworkXNoPath:
            return (None, None, None,
                    "No path found through tower network (even after bridging the nearest gap)",
                    from_ring, to_ring, None)

    plen = nx.shortest_path_length(g_conn, "from_sub", "to_sub", weight="weight")
    path_coords = []
    for node in path:
        if node == "from_sub":
            path_coords.append((from_lon, from_lat))
        elif node == "to_sub":
            path_coords.append((to_lon, to_lat))
        else:
            path_coords.append(tower_coords[node])
    return path_coords, LineString(path_coords), plen, None, from_ring, to_ring, bridge_km


# --------------------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------------------

def main():
    run_start = time.time()

    if not os.path.exists(FAILED_XLSX):
        raise FileNotFoundError(f"Failed connections file not found: {FAILED_XLSX}")
    if not os.path.exists(TOWER_SHP):
        raise FileNotFoundError(f"Tower shapefile not found: {TOWER_SHP}")

    print("=" * 80)
    print(f"ROUTING FAILED CONNECTIONS (tower-tower <= {MAX_TOWER_KM} km, "
          f"substation search 0-{MAX_SUBSTATION_SEARCH_KM} km in {SUBSTATION_RING_STEP_KM} km rings, "
          f"bridge fallback up to {BRIDGE_MAX_KM} km)")
    print("=" * 80)

    print(f"\nLoading failed connections: {FAILED_XLSX}")
    df = pd.read_excel(FAILED_XLSX)
    print(f"Loaded {len(df)} failed connections")

    if ROW_LIMIT is not None and ROW_LIMIT < len(df):
        print(f"ROW_LIMIT is set to {ROW_LIMIT} -- only processing the first {ROW_LIMIT} rows "
              f"as a test batch. Set ROW_LIMIT = None at the top of this script (or run with "
              f"UP3_ROW_LIMIT=0) to process all {len(df)} rows.")
        df = df.head(ROW_LIMIT).copy()

    print(f"\nLoading tower data: {TOWER_SHP}")
    towers_gdf = gpd.read_file(TOWER_SHP)
    print(f"Loaded {len(towers_gdf)} towers")

    if "Lat" not in towers_gdf.columns or "Long" not in towers_gdf.columns:
        print("Tower columns:", list(towers_gdf.columns))
        raise ValueError("Tower shapefile must have 'Lat' and 'Long' columns")

    required_cols = ["Fr_SS_Lat", "Fr_SS_Long", "To_SS_Lat", "To_SS_Long"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print("Excel columns:", list(df.columns))
        raise ValueError(f"Failed connections file must have these columns: {missing}")

    t0 = time.time()
    print("\nBuilding tower network...")
    g, tower_coords = build_tower_network(towers_gdf)
    print(f"Network built: {g.number_of_nodes()} nodes, {g.number_of_edges()} edges ({time.time()-t0:.1f}s)")

    t0 = time.time()
    print("Building substation-lookup spatial index...")
    tower_lons = towers_gdf["Long"].values
    tower_lats = towers_gdf["Lat"].values
    tower_tree = cKDTree(np.column_stack([tower_lons, tower_lats]))
    print(f"Spatial index built ({time.time()-t0:.1f}s)")

    print(f"\nProcessing {len(df)} failed connections...")
    print("=" * 80)

    all_lines: list[dict] = []
    failed_connections: list[dict] = []
    bridged_count = 0

    for idx, row in df.iterrows():
        row_start = time.time()
        try:
            from_lat, from_lon = row["Fr_SS_Lat"], row["Fr_SS_Long"]
            to_lat, to_lon = row["To_SS_Lat"], row["To_SS_Long"]
            circuits = str(row.get("No_of_ckt", "S/C")).strip().upper()

            if pd.isna(from_lat) or pd.isna(from_lon) or pd.isna(to_lat) or pd.isna(to_lon):
                error_msg = f"Row {idx+1}: Missing coordinates"
                print(f"! {error_msg}")
                failed_connections.append({**row.to_dict(), "Error": error_msg})
                continue
            if str(from_lon) == "#N/A" or str(to_lon) == "#N/A":
                error_msg = f"Row {idx+1}: Invalid coordinate (#N/A)"
                print(f"! {error_msg}")
                failed_connections.append({**row.to_dict(), "Error": error_msg})
                continue

            from_lat, from_lon = float(from_lat), float(from_lon)
            to_lat, to_lon = float(to_lat), float(to_lon)

            from_ss = row.get("Fr_SS", "N/A")
            to_ss = row.get("To_SS", "N/A")
            print(f"\nProcessing row {idx+1}/{len(df)}: {from_ss} -> {to_ss}", flush=True)

            (path_coords, path_line, total_distance, err,
             from_ring, to_ring, bridge_km) = connect_substations_from_coords(
                from_lon, from_lat, to_lon, to_lat,
                towers_gdf, tower_tree, tower_lons, tower_lats,
                g, tower_coords,
            )
            row_elapsed = time.time() - row_start

            if err:
                print(f"  Failed: {err}  ({row_elapsed:.1f}s)", flush=True)
                failed_connections.append({**row.to_dict(), "Error": err})
                continue

            bridge_note = f", BRIDGED {bridge_km:.2f}km gap" if bridge_km else ""
            print(f"  Path: {total_distance:.2f} km, {len(path_coords)} points "
                  f"(from-substation tower found within {from_ring:.0f}km, "
                  f"to-substation within {to_ring:.0f}km{bridge_note})  ({row_elapsed:.1f}s)", flush=True)
            if bridge_km:
                bridged_count += 1

            straight_km = haversine_distance(from_lon, from_lat, to_lon, to_lat)
            if straight_km > 0 and total_distance > straight_km * 3 and total_distance > 20:
                print(f"  Note: routed distance is {total_distance/straight_km:.1f}x the straight-line "
                      f"distance ({straight_km:.1f} km) -- worth a manual look", flush=True)

            bridge_val = round(bridge_km, 3) if bridge_km else 0.0

            if circuits == "S/C":
                line_data = row.to_dict()
                line_data["distance_km"] = total_distance
                line_data["num_towers"] = len(path_coords) - 2
                line_data["bridge_km"] = bridge_val
                all_lines.append({"geometry": path_line, **line_data})
            elif circuits == "D/C":
                line1_data = row.to_dict()
                line1_data["distance_km"] = total_distance
                line1_data["num_towers"] = len(path_coords) - 2
                line1_data["bridge_km"] = bridge_val
                line1_data["circuit"] = "Line1"
                all_lines.append({"geometry": path_line, **line1_data})
                try:
                    offset_line = create_parallel_line(path_line, offset_meters=10)
                    line2_data = row.to_dict()
                    line2_data["distance_km"] = total_distance
                    line2_data["num_towers"] = len(path_coords) - 2
                    line2_data["bridge_km"] = bridge_val
                    line2_data["circuit"] = "Line2"
                    all_lines.append({"geometry": offset_line, **line2_data})
                    print("  Parallel line (10 m offset)")
                except Exception as e:
                    print(f"  Warning: Could not create parallel line: {e}")
                    line2_data = row.to_dict()
                    line2_data["distance_km"] = total_distance
                    line2_data["num_towers"] = len(path_coords) - 2
                    line2_data["bridge_km"] = bridge_val
                    line2_data["circuit"] = "Line2"
                    all_lines.append({"geometry": path_line, **line2_data})
            else:
                print(f"  Warning: Unknown circuit type {circuits!r}, using single line")
                line_data = row.to_dict()
                line_data["distance_km"] = total_distance
                line_data["num_towers"] = len(path_coords) - 2
                line_data["bridge_km"] = bridge_val
                all_lines.append({"geometry": path_line, **line_data})

        except Exception as e:
            row_elapsed = time.time() - row_start
            error_msg = f"Row {idx+1}: {e!s}"
            print(f"  Error: {error_msg}  ({row_elapsed:.1f}s)", flush=True)
            failed_connections.append({**row.to_dict(), "Error": error_msg})
            continue

    print(f"\n{'='*80}")
    print("Creating output shapefile...")
    print(f"Successfully processed: {len(all_lines)} lines")
    print(f"  of which auto-bridged a data gap: {bridged_count}")
    print(f"Failed connections: {len(failed_connections)}")

    if all_lines:
        output_gdf = gpd.GeoDataFrame(all_lines, crs=towers_gdf.crs)
        base = OUTPUT_SHP.replace(".shp", "")
        os.makedirs(os.path.dirname(OUTPUT_SHP), exist_ok=True)
        for ext in (".shp", ".shx", ".dbf", ".prj", ".cpg", ".sbn", ".sbx"):
            p = base + ext
            if os.path.exists(p):
                try:
                    os.remove(p)
                except OSError:
                    pass
        output_gdf.to_file(OUTPUT_SHP)
        print(f"\nLines saved: {OUTPUT_SHP} ({len(output_gdf)} features)")
        print("  Check the 'bridge_km' column -- any value > 0 means that line contains "
              "one auto-bridged stretch across a data gap, worth a manual look.")
    else:
        print("\nNo lines were created.")

    if failed_connections:
        failed_df = pd.DataFrame(failed_connections)
        os.makedirs(os.path.dirname(REMAINING_FAILED_XLSX), exist_ok=True)
        failed_df.to_excel(REMAINING_FAILED_XLSX, index=False)
        print(f"\nFailed connections: {REMAINING_FAILED_XLSX} ({len(failed_connections)} rows)")
    else:
        print("\nAll connections routed.")

    total_elapsed = time.time() - run_start
    print(f"\n{'='*80}")
    print(f"Done. Total time: {total_elapsed/60:.1f} minutes")
    print(f"{'='*80}")


def run_with_logging():
    project_dir = _script_dir()
    os.chdir(project_dir)
    log_dir = project_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = log_dir / f"route_v3_{timestamp}.log"

    real_stdout, real_stderr = sys.stdout, sys.stderr
    with log_path.open("w", encoding="utf-8") as log_file:
        sys.stdout = Tee(real_stdout, log_file)
        sys.stderr = Tee(real_stderr, log_file)
        try:
            print(f"Logging to: {log_path}\n")
            try:
                main()
            except Exception:
                import traceback
                traceback.print_exc()
                print(f"\nRun failed. Full log saved to: {log_path}")
                raise
            else:
                print(f"\nFull log saved to: {log_path}")
        finally:
            sys.stdout = real_stdout
            sys.stderr = real_stderr


if __name__ == "__main__":
    run_with_logging()
