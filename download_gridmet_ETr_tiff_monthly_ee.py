import ee
# import sys
import pandas as pd
import os
import math
import logging
import sys

ee.Initialize()




# Start Date, End Date (Exclusive), and List for File Naming
# Year Filter
year = 2017

# Constants
pi = math.pi

# Read in/calculate ancillary layers
elev = ee.Image("projects/climate-engine/gridmet/elevation")
lat = ee.Image.pixelLonLat().select('latitude').multiply(pi / 180)
# lon = ee.Image.pixelLonLat().select('longitude').multiply(pi / 180)
pair = elev.expression('101.3 * pow((293 - 0.0065 * b()) / 293, 5.26)')
# pair = elev.expression('101.3 * pow((293 - 0.0065 * b()) / 293, 5.26)')

# GRIDMET Daily ETr
# rhmin = gridmet_image.select(["rmin"]).multiply(0.01)  # % to decimal
# rhmax = gridmet_image.select(["rmax"]).multiply(0.01)  # % to decimal
# Vapor pressure from RHmax and RHmin (Eqn 11)
# ea = es_tmin.multiply(rhmax).add(es_tmax.multiply(rhmin)).multiply(0.5)
# Vapor pressure from specific humidity (Eqn )
# To match standardized form, ea is calculated from elevation based pair
def gridmet_etr_func (gridmet_image):
    scene_date = ee.Algorithms.Date(gridmet_image.get("system:time_start"))
    doy = ee.Number(scene_date.getRelative('day', 'year')).add(1).double()
    tmin = gridmet_image.select(["tmmn"]).subtract(273.15)
    tmax = gridmet_image.select(["tmmx"]).subtract(273.15)
    q = gridmet_image.select(["sph"])
    rs = gridmet_image.select(["srad"]).multiply(0.0864)
    uz = gridmet_image.select(["vs"])
    zw = 10.0
    ea = pair.expression('q * pair / (0.622 + 0.378 * q)', {'pair': pair, 'q': q})
    return daily_refet_func(doy, tmin, tmax, ea, rs, uz, zw, 1600, 0.38) \
    .select([0], ['ETr']) \
    .copyProperties(gridmet_image, ['system:index', 'system:time_start', 'system:time_end'])

# Saturated Vapor Pressure in kPa with temperature in C
def vapor_pressure_func (t):
    return t.expression('0.6108 * exp(17.27 * b() / (b() + 237.3))')

# Daily Reference ET
def daily_refet_func(doy, tmin, tmax, ea, rs, uz, zw, cn, cd):
    # Calculations
    tmean = tmin.add(tmax).multiply(0.5)  # C
    # To match standardized form, psy is calculated from elevation based pair
    psy = pair.multiply(0.000665)
    es_tmax = vapor_pressure_func(tmax)  # C
    es_tmin = vapor_pressure_func(tmin)  # C
    es_tmean = vapor_pressure_func(tmean)
    es_slope = es_tmean.expression(
        '4098 * es / (pow((t + 237.3), 2))', {'es': es_tmean, 't': tmean})
    es = es_tmin.add(es_tmax).multiply(0.5)
    # Extraterrestrial radiation (Eqn 24, 27, 23, 21)
    delta = ee.Image.constant(
        doy.multiply(2 * pi / 365).subtract(1.39).sin().multiply(0.409))
    # delta = ee.Image.constant(
    #  doy.multiply(2 * pi / 365).subtract(1.39435).sin().multiply(0.40928))
    omegas = lat.expression(
        'acos(-tan(lat) * tan(delta))', {'lat': lat, 'delta': delta})
    theta = omegas.expression(
        'omegas * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(b())',
        {'omegas': omegas, 'lat': lat, 'delta': delta})
    dr = ee.Image.constant(
        doy.multiply(2 * pi / 365).cos().multiply(0.033).add(1))
    ra = theta.expression(
        '(24 / pi) * gsc * dr * theta',
        {'pi': pi, 'gsc': 4.92, 'dr': dr, 'theta': theta})

    # Simplified clear sky solar formulation (Eqn 19)
    # rso = elev.expression(
    #  '(0.75 + 2E-5 * elev) * ra', {'elev':elev, 'ra':ra})
    # addToMap(rso, {}, 'rso')

    # This is the full clear sky solar formulation
    # sin of the angle of the sun above the horizon (D.5 and Eqn 62)
    sin_beta_24 = lat.expression(
        'sin(0.85 + 0.3 * lat * sin(delta) - 0.42 * lat ** 2)',
        {'lat': lat, 'delta': ee.Image.constant(doy.multiply(2 * pi / 365).subtract(1.39435))})
    # Precipitable water (Eqn D.3)
    w = pair.expression('0.14 * ea * pair + 2.1', {'pair': pair, 'ea': ea})
    # Clearness index for direct beam radiation (Eqn D.2)
    # Limit sin_beta >= 0.01 so that KB does not go undefined
    kb = pair.expression(
        '0.98 * exp((-0.00146 * pair) / (kt * sin_beta) - 0.075 * pow((w / sin_beta), 0.4))',
        {'pair': pair, 'kt': 1.0, 'sin_beta': sin_beta_24.max(0.01), 'w': w})
    # Transmissivity index for diffuse radiation (Eqn D.4)
    kd = kb.multiply(-0.36).add(0.35) \
    .min(kb.multiply(0.82).add(0.18))
    # kd = kb.multiply(-0.36).add(0.35)
    #  .where(kb.lt(0.15), kb.multiply(0.82).add(0.18))
    # (Eqn D.1)
    rso = ra.multiply(kb.add(kd))

    # Cloudiness fraction (Eqn 18)
    fcd = rs.divide(rso).clamp(0.3, 1).multiply(1.35).subtract(0.35)

    # Net long-wave radiation (Eqn 17)
    rnl = ea.expression(
        '4.901E-9 * fcd * (0.34 - 0.14 * sqrt(ea)) * ' +
        '(pow(tmax_k, 4) + pow(tmin_k, 4)) / 2',
        {'ea': ea, 'fcd': fcd, 'tmax_k': tmax.add(273.15),
         'tmin_k': tmin.add(273.15)})

    # Net radiation (Eqns 15 and 16)
    rn = rs.multiply(0.77).subtract(rnl)

    # Wind speed (Eqn 33)
    u2 = uz.expression('b() * 4.87 / log(67.8 * zw - 5.42)', {'zw': zw})

    # Daily reference ET (Eqn 1)
    return tmin.expression(
        ('(0.408 * slope * rn + (psy * cn * u2 * (es - ea) / (t + 273))) / ' +
         '(slope + psy * (cd * u2 + 1))'),
        {'slope': es_slope, 'rn': rn, 'psy': psy, 'cn': cn,
         't': tmean, 'u2': u2, 'es': es, 'ea': ea, 'cd': cd})



    # Clip Boundary
# HUC2 = ee.FeatureCollection('USGS/WBD/2017/HUC02')
# watershed = HUC2.filterMetadata('name', 'equals', 'Upper Colorado Region').geometry().bounds()
# Transform taken from GRIDMET ee.image
# var gridmet_transform = JSON.stringify(gridmet.projection().getInfo()['transform']);
gridmet_transform = [0.041666666666666664,0,-124.78749996666667,0,-0.041666666666666664,49.42083333333334]

for month in range(1, 13):
    gridmet_coll = ee.ImageCollection('IDAHO_EPSCOR/GRIDMET') \
        .filter(ee.Filter.calendarRange(year, year, 'year')) \
        .filter(ee.Filter.calendarRange(month, month, 'month')) \

    gridmet_etr = gridmet_coll.map(gridmet_etr_func)

    gridmet_image = ee.Image(gridmet_etr.sum())

    # Map month_stat_fnc over month list. Create collection from images.
    ee.batch.Export.image.toCloudStorage(
        image = gridmet_image,
        bucket = 'ee_exports',
        fileNamePrefix = '{}/{}_Sum_{}_{:02d}'.format('ETr', 'ETr', year, month),
        dimensions ='1386x585',
        crs = 'EPSG:4326',
        crsTransform = gridmet_transform).start()


