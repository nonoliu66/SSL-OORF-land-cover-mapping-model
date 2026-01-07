/************************************************************
 * Integrated Feature Stack for SNIC Segmentation
 *
 * Feature components:
 * 1) NDVI-based phenological features (Table S2)
 * 2) Sentinel-2 spectral indices (Table S1)
 * 3) GLCM texture features (Table S3)
 * 4) Sentinel-2 median reflectance bands
 * 5) Sentinel-1 VV / VH mean backscatter
 *
 * Output:
 * - SNIC clusters with red boundaries (transparent interiors)
 *
 * Required assets:
 * - table   : study area
 * - table2 : reference points (optional)
 ************************************************************/


// ==========================================================
// 0. Study area
// ==========================================================
var roi = table.geometry().bounds();
Map.centerObject(roi, 11);


// ==========================================================
// 1. NDVI-based phenological features (Table S2)
// ==========================================================
var year = '2024';
var startDate = ee.Date(year + '-01-01');
var endDate   = ee.Date(year + '-12-31');

// Monthly NDVI time series
var s2_ndvi = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(function(img) {
    var ndvi = img.normalizedDifference(['B8', 'B4'])
                  .rename('NDVI');
    return ndvi.copyProperties(img, ['system:time_start']);
  });

var months = ee.List.sequence(1, 12);

// Monthly composites
var ndviMonthly = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var start = startDate.advance(ee.Number(m).subtract(1), 'month');
    var end   = start.advance(1, 'month');
    return s2_ndvi.filterDate(start, end)
                  .mean()
                  .set('month', m)
                  .set('system:time_start', start.millis());
  })
);

var ndviList = ndviMonthly.toList(12);

// NDVI maximum and threshold (50% of NDVImax)
var NDVImax = ndviMonthly.max().rename('NDVImax');
var NDVIthr = NDVImax.multiply(0.5).rename('NDVIthr');

// Growing season timing
var monthMaskCol = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    return ee.Image.constant(m)
      .updateMask(nd.gt(NDVIthr))
      .toFloat();
  })
);

var Month_start_val = monthMaskCol.min().unmask(0).rename('Month_start_val');
var Month_end_val   = monthMaskCol.max().unmask(0).rename('Month_end_val');

// NDVI at start and end of season
var Start_val = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    return nd.updateMask(Month_start_val.eq(m));
  })
).mosaic().unmask(0).rename('Start_val');

var End_val = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    return nd.updateMask(Month_end_val.eq(m));
  })
).mosaic().unmask(0).rename('End_val');

// Base value and peak value
var Base_val = Start_val.add(End_val).divide(2).rename('Base_val');

var growingSeasonCol = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    var cond = ee.Image.constant(m)
      .gte(Month_start_val)
      .and(ee.Image.constant(m).lte(Month_end_val));
    return nd.updateMask(cond);
  })
);

var Peak_val = growingSeasonCol.max().rename('Peak_val');
var Ampl     = Peak_val.subtract(Base_val).rename('Ampl');

// Month of peak NDVI
var tol = 1e-6;
var Month_peak_val = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    return ee.Image.constant(m)
      .updateMask(nd.subtract(Peak_val).abs().lte(tol))
      .toFloat();
  })
).min().unmask(0).rename('Month_peak_val');

// Seasonal cumulative NDVI
var L_cumulative = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    return nd.updateMask(
      ee.Image.constant(m).gte(Month_start_val)
        .and(ee.Image.constant(m).lte(Month_peak_val))
    );
  })
).sum().unmask(0).rename('L_cumulative');

var R_cumulative = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    return nd.updateMask(
      ee.Image.constant(m).gte(Month_peak_val)
        .and(ee.Image.constant(m).lte(Month_end_val))
    );
  })
).sum().unmask(0).rename('R_cumulative');

// NDVI change rates
var diffCol = ee.ImageCollection.fromImages(
  ee.List.sequence(1, 11).map(function(i) {
    var before = ee.Image(ndviList.get(ee.Number(i).subtract(1)));
    var after  = ee.Image(ndviList.get(i));
    return after.subtract(before);
  })
);

var Max_increase = diffCol.max().rename('Max_increase');
var Min_decrease = diffCol.min().multiply(-1).rename('Min_decrease');

// Phenological feature stack
var phenology = NDVImax.addBands([
  NDVIthr,
  Base_val,
  Start_val,
  End_val,
  Peak_val,
  Ampl,
  Month_start_val,
  Month_peak_val,
  Month_end_val,
  L_cumulative,
  R_cumulative,
  Max_increase,
  Min_decrease
]).clip(roi).toFloat();


// ==========================================================
// 2. Sentinel-2 spectral indices (Table S1)
// ==========================================================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(function(img) {
    return img.divide(10000)
      .select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'])
      .copyProperties(img, ['system:time_start']);
  });

var s2_year = s2.median().clip(roi);

// Spectral indices
var NDVI  = s2_year.normalizedDifference(['B8','B4']).rename('NDVI');
var NDWI  = s2_year.normalizedDifference(['B3','B8']).rename('NDWI');
var GCVI  = s2_year.expression('B8/B3 - 1', {'B8':s2_year.select('B8'),'B3':s2_year.select('B3')}).rename('GCVI');
var EVI   = s2_year.expression('2.5*(B8-B4)/(B8+6*B4-7.5*B2+1)', {'B8':s2_year.select('B8'),'B4':s2_year.select('B4'),'B2':s2_year.select('B2')}).rename('EVI');
var SAVI  = s2_year.expression('((B8-B4)/(B8+B4+0.5))*1.5', {'B8':s2_year.select('B8'),'B4':s2_year.select('B4')}).rename('SAVI');
var RVI   = s2_year.expression('B8/B4', {'B8':s2_year.select('B8'),'B4':s2_year.select('B4')}).rename('RVI');
var WI    = s2_year.expression('(B3+2.5*B11)/(1.5*B8+2*B4+0.5*B11)', {'B3':s2_year.select('B3'),'B11':s2_year.select('B11'),'B8':s2_year.select('B8'),'B4':s2_year.select('B4')}).rename('WI');
var MSI   = s2_year.expression('B11/B8', {'B11':s2_year.select('B11'),'B8':s2_year.select('B8')}).rename('MSI');
var NDBI  = s2_year.normalizedDifference(['B11','B8']).rename('NDBI');
var NDMI  = s2_year.normalizedDifference(['B8','B11']).rename('NDMI');
var MNDWI = s2_year.normalizedDifference(['B3','B11']).rename('MNDWI');
var DVI   = s2_year.expression('B8-B4', {'B8':s2_year.select('B8'),'B4':s2_year.select('B4')}).rename('DVI');

var spectral = NDVI.addBands([
  NDWI, GCVI, EVI, SAVI, RVI, WI, MSI, NDBI, NDMI, MNDWI, DVI
]);


// ==========================================================
// 3. GLCM texture features (Table S3)
// ==========================================================
var gray = s2_year.expression(
  '0.3*B8 + 0.59*B4 + 0.11*B3',
  {'B8':s2_year.select('B8'),'B4':s2_year.select('B4'),'B3':s2_year.select('B3')}
).rename('Gray');

var glcm = gray.multiply(100).toInt().glcmTexture({size:3});

var texture = glcm.select(
  ['Gray_contrast','Gray_asm','Gray_corr','Gray_var',
   'Gray_idm','Gray_diss','Gray_sent','Gray_ent'],
  ['gray_contrast','gray_asm','gray_corr','gray_var',
   'gray_idm','gray_diss','gray_sent','gray_ent']
);


// ==========================================================
// 4. Sentinel-1 SAR backscatter features
// ==========================================================
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));

var VV_mean = s1.select('VV').mean().rename('VV');
var VH_mean = s1.select('VH').mean().rename('VH');


// ==========================================================
// 5. Final feature stack for SNIC
// ==========================================================
var allFeatures = phenology
  .addBands(spectral)
  .addBands(texture)
  .addBands(s2_year)
  .addBands([VV_mean, VH_mean])
  .clip(roi)
  .toFloat();

print('Feature list:', allFeatures.bandNames());


// ==========================================================
// 6. SNIC segmentation and red boundary visualization
// ==========================================================
function runSNIC(image, seedSpacing, size) {
  var seeds = ee.Algorithms.Image.Segmentation.seedGrid(seedSpacing);
  return ee.Algorithms.Image.Segmentation.SNIC({
    image: image,
    size: size,
    compactness: 5,
    connectivity: 8,
    neighborhoodSize: 2 * size + 1,
    seeds: seeds
  }).select('clusters');
}

function extractRedEdges(clusters) {
  var minNb = clusters.reduceNeighborhood({
    reducer: ee.Reducer.min(),
    kernel: ee.Kernel.square(1)
  });
  return minNb.neq(clusters)
    .selfMask()
    .visualize({palette:['FF0000'], opacity:1});
}

// Example segmentation
var imageForSNIC = allFeatures.select([
  'Peak_val','Base_val','NDVI','RVI','MNDWI',
  'Start_val','End_val','B2','B3','B4','B11','B12',
  'VV','VH'
]);

var seg40 = runSNIC(imageForSNIC, 40, 40);
Map.addLayer(extractRedEdges(seg40), {}, 'SNIC edges (40)');

// Multi-scale segmentation
var scales = [20, 40, 60, 80, 100];
scales.forEach(function(s) {
  var seg = runSNIC(allFeatures, s, s);
  Map.addLayer(extractRedEdges(seg), {}, 'edges_' + s);
});

// Export SNIC results
scales.forEach(function(s) {
  var seg = runSNIC(allFeatures, s, s);
  Export.image.toDrive({
    image: seg,
    description: 'SNIC_seg_' + s,
    folder: 'GEE_Export',
    scale: 10,
    region: roi,
    maxPixels: 1e13
  });
});
