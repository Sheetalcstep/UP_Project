# UP_Project

Scripts for routing Uttar Pradesh power transmission lines through the power-tower
network (from substation coordinates, via existing towers, to another substation).

## Index

1. [Overview](#overview)
2. [Getting the code](#getting-the-code)
3. [Files in this repo](#files-in-this-repo)
4. [Transmission_line_UP.py](#transmission_line_uppy)
   - [Purpose](#purpose)
   - [Distances used](#distances-used)
   - [Inputs](#inputs)
   - [Outputs](#outputs)
   - [How to run](#how-to-run)
5. [route_failed_connections_v3.py](#route_failed_connections_v3py)
   - [Purpose](#purpose-1)
   - [What's different from Transmission_line_UP.py](#whats-different-from-transmission_line_uppy)
   - [Distances used](#distances-used-1)
   - [Bridge fallback for data gaps](#bridge-fallback-for-data-gaps)
   - [Inputs](#inputs-1)
   - [Outputs](#outputs-1)
   - [How to run](#how-to-run-1)
   - [Config / environment variables](#config--environment-variables)
6. [Requirements](#requirements)
7. [Known limitations](#known-limitations)

## Overview

Both scripts build a graph of power towers (from a point shapefile), connect each
substation to nearby towers, and use a shortest-path search through the tower
network to draw a line between two substations -- instead of a straight
substation-to-substation line that ignores the real tower trail.

## Getting the code

Clone the repo and move into it:

```
git clone https://github.com/Sheetalcstep/UP_Project.git
cd UP_Project
```

## Files in this repo

| File | Purpose |
|---|---|
| `Transmission_line_UP.py` | Original, general-purpose routing script (all connections, from a CSV). Includes water-crossing bridge logic. |
| `route_failed_connections_v3.py` | Routing script specifically for the *failed* connections list, with a wider substation search radius and an automatic gap-bridging fallback. |

## Transmission_line_UP.py

### Purpose

Batch-routes substation-to-substation connections listed in a CSV, using the
power-tower network as the path, and avoiding water crossings where possible.

### Distances used

- Tower-to-tower: 0.5 km (500 m) on land by default.
- Tower-to-tower: up to 1.2 km allowed *only* if the straight segment between
  the two towers crosses a water feature (`water_bodies.shp`) -- a wider span
  just to bridge rivers/lakes where towers are naturally spaced further apart.
- Substation-to-tower: 1.0 km -- a substation connects to every tower within
  1 km (not just the closest one); the shortest-path search picks whichever
  entry point gives the best route.

### Inputs

- `excel/UP_power_map_line_length_csv_SC_DC.csv` -- connections to route
  (columns: `Fr_SS_Lat`, `Fr_SS_Long`, `To_SS_Lat`, `To_SS_Long`, `No_of_ckt`, ...)
- `power_substation_UP_other_2.shp` -- tower points (`Lat`, `Long` columns)
- `water_bodies.shp` -- water polygons, used for the water-crossing bridge rule

### Outputs

- `Shapefile/UP_transmission_lines.shp` -- routed lines
- `excel/failed_central_UP.xlsx` -- connections that could not be routed

### How to run

```
python Transmission_line_UP.py
```

File paths can be overridden with the `UP_CSV`, `UP_TOWER_SHP`, `UP_WATER_SHP`,
`UP_OUT_SHP`, and `UP_FAILED_XLSX` environment variables.

## route_failed_connections_v3.py

### Purpose

Routes specifically the connections that `Transmission_line_UP.py` (or an
earlier pass) could not resolve -- most of which failed either because no
tower was within 1 km of a substation, or because of a genuine gap in the
tower network's connectivity.

### What's different from Transmission_line_UP.py

- No water-crossing exception -- tower-to-tower is a flat distance, no
  `water_bodies.shp` dependency.
- Substation-to-tower search expands outward in rings (1 km, 2 km, 3 km, ...)
  instead of a single fixed radius, since these rows specifically failed due
  to substations being further than 1 km from the nearest tower.
- Adds a bridge fallback for cases where a substation's nearby towers form
  their own small, disconnected island in the network (see below) -- these
  rows still failed under the original approach even with a wider search
  radius, because the *tower network itself* had a gap, not the substation
  search.

### Distances used

- Tower-to-tower: 1 km flat (`MAX_TOWER_KM`).
- Substation-to-tower: rings of 1 km (`SUBSTATION_RING_STEP_KM`), expanding up
  to 30 km (`MAX_SUBSTATION_SEARCH_KM`). Stops at the first ring that contains
  at least one tower and connects to every tower in that ring.
- Bridge fallback cap: 8 km (`BRIDGE_MAX_KM`) -- see below.

### Bridge fallback for data gaps

Some rows fail not because a substation has no nearby tower, but because the
towers it finds sit in their own small island, cut off from the rest of the
tower network by a real gap in the source data (missing towers along that
stretch). Widening the 1 km tower-to-tower rule for everyone risks reconnecting
things that shouldn't be connected elsewhere, so instead:

1. Try the normal path first.
2. If that fails, find the two tower network "islands" (connected components)
   the from/to substations actually landed in, and measure the real gap
   between them (closest tower-to-tower distance across the two islands).
3. If that gap is <= `BRIDGE_MAX_KM` (default 8 km), add one straight bridge
   segment across just that gap, reconnect, and re-route. The output row is
   tagged with a `bridge_km` value so you can filter for exactly which lines
   contain an auto-bridged stretch and check them on the map.
4. If the gap is bigger than `BRIDGE_MAX_KM`, the row is left as failed for
   manual review, and the error message reports how big the real gap is.

### Inputs

- `excel/failed_UP_3_remaining.xlsx` -- the still-failing connections list
  (same columns as above)
- `power_substation_UP_other_2.shp` -- tower points

### Outputs

- `Shapefile/failed_lines_UP_bridged.shp` -- routed lines, including a
  `bridge_km` column (0 if no bridge was needed)
- `excel/failed_UP_bridged_remaining.xlsx` -- connections still unresolved
  even after the bridge fallback
- `logs/route_v3_<timestamp>.log` -- full console output of the run

### How to run

```
python route_failed_connections_v3.py
```

`ROW_LIMIT` near the top of the script controls how many rows get processed
(default: 10, for a quick test batch). Set it to `None` in the file, or run
with `UP3_ROW_LIMIT=0`, to process the full file.

### Config / environment variables

| Variable | Default | Purpose |
|---|---|---|
| `UP3_FAILED_XLSX` | `excel/failed_UP_3_remaining.xlsx` | Input connections list |
| `UP3_TOWER_SHP` | `power_substation_UP_other_2.shp` | Tower points |
| `UP3_OUT_SHP` | `Shapefile/failed_lines_UP_bridged.shp` | Output lines |
| `UP3_REMAINING_XLSX` | `excel/failed_UP_bridged_remaining.xlsx` | Still-failed rows |
| `UP3_ROW_LIMIT` | `10` (script default) | Rows to process; `0` = all rows |

## Requirements

- Python 3.10+
- `geopandas`, `pandas`, `networkx`, `shapely`, `scipy`, `numpy`, `pyproj`

## Known limitations

- Tower data comes from OpenStreetMap points with no attribute indicating
  which physical line/circuit a tower belongs to -- routing is purely
  distance-based, so results are only as good as the tower coverage in the
  source data.
- The bridge fallback patches genuine data gaps with a single straight
  segment; lines with `bridge_km > 0` should be spot-checked on a map before
  being treated as final.
