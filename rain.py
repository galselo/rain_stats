# WARNING: this code is for reference only, use at your own risk
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import datetime
from scipy.interpolate import interp2d

# download data at
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/insitu-gridded-observations-europe?tab=overview


# station names with lat, lon, mm/day
stations = {"Rheinbach-Todenfeld": [50.58, 6.948, 158.],
            "Cologne-Stammheim": [50.989, 6.978, 154.],
            "Klein-Altendorf": [50.623, 6.992, 147.],
            "Kall-Sistig": [50.501, 6.526, 145.]}

# extract data for later use
station_names = sorted(stations.keys())
station_lats = [stations[k][0] for k in station_names]
station_lons = [stations[k][1] for k in station_names]
station_mms = [stations[k][2] for k in station_names]

print(station_lats)

# read precipitation dataset
f = Dataset("rr_ens_mean_0.1deg_reg_v23.1e.nc")

# get number of time-steps (days)
ntimes = f["time"].shape[0]

# extract latitued and longitude for later use
lat = np.array(f["latitude"])
lon = np.array(f["longitude"])

# store year days for specific range (e.g. JAS)
days = []
for y in range(1950, 1951):
  d1 = datetime.date(1950, 1, 1)  # initial date in the dataset

  d2 = datetime.date(y, 7, 1)  # starting date in the year
  d3 = datetime.date(y, 10, 1)  # ending date in the year

  # compute the range from the beginning of the dataset
  days += range((d2-d1).days, (d3-d1).days+1)


# latitude and longitude range
latmin = 49.5
latmax = 51.5
lonmin = 5.
lonmax = 9.

# number of days to sum, e.g. ndays=2 is mm/48h
ndays = 1

# find size of the grid to initialize stats arrays
rr = f["rr"][0, (latmin < lat) & (lat < latmax), (lonmin < lon) & (lon < lonmax)]
xx, yy = rr.shape
st = np.zeros((len(days), xx, yy))
mx = np.zeros_like(rr)
rrtot = f["rr"][0, ...]
mxtot = np.zeros_like(rrtot)


# loop on selected days
for i, d in enumerate(tqdm(days)):
  # sum over days (d+1 means no sum, i.e. single day)
  rr = np.sum(f["rr"][d:d+ndays, (latmin < lat) & (lat < latmax), (lonmin < lon) & (lon < lonmax)], axis=0)

  # use full map for stats
  rrtot = np.sum(f["rr"][d:d+ndays, ...], axis=0)

  # add data to statistics matrix
  st[i, ...] = rr
  # find maximum
  mx = np.maximum(rr, mx)
  mxtot = np.maximum(rrtot, mxtot)

print("maximum precipitation:", np.amax(st))

# flatten data and remove negative data (i.e. no data)
r = st.flatten()
r = r[r >= 0]

# plot historgram of data
plt.hist(r, alpha=0.2, bins=30, density=True)

# add stations to hist
for mm in station_mms:
  plt.axvline(mm, color="k", ls="--")


# additional stuff
plt.yscale("log")
plt.xlabel("mm/%dh" % (ndays*24))
plt.ylabel("Normalized frequency")
plt.savefig("hist.png", dpi=180)
print("histogram saved as hist.png")


plt.clf()
# prepare lat/long grid for map plot
xv, yv = np.meshgrid(lon[(lonmin < lon) & (lon < lonmax)], lat[(latmin < lat) & (lat < latmax)], indexing="ij")

plt.pcolor(xv, yv, mx.T, cmap="cubehelix")
plt.colorbar()
plt.scatter(station_lons, station_lats, marker="x", color="w", s=100)

for llt, lln, mm in zip(station_lats, station_lons, station_names):
  plt.text(lln+0.1, llt, mm, color="w", rotation=0, ha="left", va="center")

plt.ylabel("Decimal latitude")
plt.xlabel("Decimal longitude")
plt.title("Max precipitation 1950-2020, mm/%dh" % (ndays*24))
plt.savefig("map.png", dpi=180)
print("map saved as map.png")


# interpolate maximum to find station max
f = interp2d(lat[(latmin < lat) & (lat < latmax)], lon[(lonmin < lon) & (lon < lonmax)], mx.T)
print("Max rain station interpolated:")
for k, v in stations.items():
  print(k, f(v[0], v[1]))



