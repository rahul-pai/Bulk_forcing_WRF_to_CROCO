# Bulk_forcing_WRF_to_CROCO
Python interface that converts WRF Lambert-grid 3-hourly output into the ERA5 intermediate format read by croco_tools (make_ERA5.m). It rotates winds to east/north, takes exact 3-hour means of radiation and precipitation (handling WRF’s running-total overflows and resets), and regrids to a regular lon/lat grid. STATUS: Under testing
