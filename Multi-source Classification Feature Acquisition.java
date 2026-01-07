/************************************************************
 * Sentinel-based Feature Extraction Framework (1-Year)
 * 
 * Components:
 * 1) NDVI-based phenological features (Table S2)
 * 2) Spectral indices from Sentinel-2 (Table S1)
 * 3) Texture features based on GLCM (Table S3)
 * 4) Multi-source feature integration (Sentinel-1 + Sentinel-2)
 *
 * All feature names strictly follow Tables S1–S3 in the paper
 ************************************************************/


// ==========================================================
// 1. General parameters
// ==========================================================
var year = '2023';
var startDate = ee.Date(year + '-01-01');
var endDate   = ee.Date(year + '-12-31');
var roi = table.geometry().bounds();   // Study area


// ==========================================================
// 2. Sentinel-2 NDVI time series (monthly)
// ==========================================================
var s2_ndvi = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(function(img) {
    var ndvi = img.normalizedDifference(['B8', 'B4'])
                  .rename('NDVI');
    return ndvi.copyProperties(img, ['system:time_start']);
  });

// Monthly NDVI composites
var months = ee.List.sequence(1, 12);

var ndviMonthly = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var start = ee.Date(year + '-01-01').advance(ee.Number(m).subtract(1), 'month');
    var end   = start.advance(1, 'month');
    return s2_ndvi.filterDate(start, end)
                  .mean()
                  .set('month', m)
                  .set('system:time_start', start.millis());
  })
);

var ndviList = ndviMonthly.toList(12);


// ==========================================================
// 3. Phenological timing based on NDVI threshold
//    Threshold = 0.5 × NDVImax (pixel-wise)
// ==========================================================
var NDVImax = ndviMonthly.max().rename('NDVImax');
var NDVIthr = NDVImax.multiply(0.5).rename('NDVIthr');

// Identify months exceeding the threshold
var monthMaskCol = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var img = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    return ee.Image.constant(m)
      .updateMask(img.gt(NDVIthr))
      .toFloat();
  })
);

// Start, peak, and end of the growing season
var Month_start_val = monthMaskCol.min().unmask(0).rename('Month_start_val');
var Month_end_val   = monthMaskCol.max().unmask(0).rename('Month_end_val');


// ==========================================================
// 4. NDVI values during the growing season
// ==========================================================

// NDVI value at the start of the season
var Start_val = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    return nd.updateMask(Month_start_val.eq(m));
  })
).mosaic().unmask(0).rename('Start_val');

// NDVI value at the end of the season
var End_val = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    return nd.updateMask(Month_end_val.eq(m));
  })
).mosaic().unmask(0).rename('End_val');

// Base value: average of Start_val and End_val
var Base_val = Start_val.add(End_val)
                        .divide(2)
                        .rename('Base_val');


// ==========================================================
// 5. Peak NDVI and seasonal amplitude
// ==========================================================
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


// Identify the month of peak NDVI
var tol = 1e-6;
var Month_peak_val = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    return ee.Image.constant(m)
      .updateMask(nd.subtract(Peak_val).abs().lte(tol))
      .toFloat();
  })
).min().unmask(0).rename('Month_peak_val');


// ==========================================================
// 6. Seasonal cumulative NDVI
// ==========================================================

// Left cumulative: start → peak
var L_cumulative = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    var cond = ee.Image.constant(m)
      .gte(Month_start_val)
      .and(ee.Image.constant(m).lte(Month_peak_val));
    return nd.updateMask(cond);
  })
).sum().unmask(0).rename('L_cumulative');

// Right cumulative: peak → end
var R_cumulative = ee.ImageCollection.fromImages(
  months.map(function(m) {
    var nd = ee.Image(ndviList.get(ee.Number(m).subtract(1)));
    var cond = ee.Image.constant(m)
      .gte(Month_peak_val)
      .and(ee.Image.constant(m).lte(Month_end_val));
    return nd.updateMask(cond);
  })
).sum().unmask(0).rename('R_cumulative');


// ==========================================================
// 7. NDVI change rates during the season
// ==========================================================
var gsList = ee.List.sequence(1, 11).map(function(i) {
  var before = ee.Image(ndviList.get(ee.Number(i).subtract(1)));
  var after  = ee.Image(ndviList.get(i));
  return after.subtract(before);
});

var diffCol = ee.ImageCollection.fromImages(gsList);

var Max_increase = diffCol.max().rename('Max_increase');
var Min_decrease = diffCol.min().multiply(-1).rename('Min_decrease');


// ==========================================================
// 8. Phenological feature stack (Table S2)
// ==========================================================
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
// 9. Spectral indices from Sentinel-2 (Table S1)
// ==========================================================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(function(img) {
    return img.divide(10000)
      .select(['B2','B3','B4','B8','B11'])
      .copyProperties(img, ['system:time_start']);
  });

var s2mean = s2.median().clip(roi);

// Spectral indices
var NDVI  = s2mean.normalizedDifference(['B8','B4']).rename('NDVI');
var NDWI  = s2mean.normalizedDifference(['B3','B8']).rename('NDWI');
var GCVI  = s2mean.expression('B8 / B3 - 1', {'B8':s2mean.select('B8'),'B3':s2mean.select('B3')}).rename('GCVI');
var EVI   = s2mean.expression('2.5*(B8-B4)/(B8+6*B4-7.5*B2+1)',{'B8':s2mean.select('B8'),'B4':s2mean.select('B4'),'B2':s2mean.select('B2')}).rename('EVI');
var SAVI  = s2mean.expression('((B8-B4)/(B8+B4+0.5))*1.5',{'B8':s2mean.select('B8'),'B4':s2mean.select('B4')}).rename('SAVI');
var RVI   = s2mean.expression('B8/B4',{'B8':s2mean.select('B8'),'B4':s2mean.select('B4')}).rename('RVI');
var WI    = s2mean.expression('(B3+2.5*B11)/(1.5*B8+2*B4+0.5*B11)',{'B3':s2mean.select('B3'),'B8':s2mean.select('B8'),'B4':s2mean.select('B4'),'B11':s2mean.select('B11')}).rename('WI');
var MSI   = s2mean.expression('B11/B8',{'B11':s2mean.select('B11'),'B8':s2mean.select('B8')}).rename('MSI');
var NDBI  = s2mean.normalizedDifference(['B11','B8']).rename('NDBI');
var NDMI  = s2mean.normalizedDifference(['B8','B11']).rename('NDMI');
var MNDWI = s2mean.normalizedDifference(['B3','B11']).rename('MNDWI');
var DVI   = s2mean.expression('B8-B4',{'B8':s2mean.select('B8'),'B4':s2mean.select('B4')}).rename('DVI');

var spectral = NDVI.addBands([
  NDWI, GCVI, EVI, SAVI, RVI, WI, MSI, NDBI, NDMI, MNDWI, DVI
]);


// ==========================================================
// 10. Texture features based on GLCM (Table S3)
// ==========================================================
var gray = s2mean.expression(
  '0.3*B8 + 0.59*B4 + 0.11*B3',
  {'B8':s2mean.select('B8'),'B4':s2mean.select('B4'),'B3':s2mean.select('B3')}
).rename('Gray');

var glcm = gray.multiply(100).toInt().glcmTexture({size:3});

var texture = glcm.select(
  ['Gray_contrast','Gray_asm','Gray_corr','Gray_var',
   'Gray_idm','Gray_diss','Gray_sent','Gray_ent'],
  ['gray_contrast','gray_asm','gray_corr','gray_var',
   'gray_idm','gray_diss','gray_sent','gray_ent']
);


// ==========================================================
// 11. Sentinel-1 SAR features
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
// 12. Final feature stack
// ==========================================================
var allFeatures = phenology
  .addBands(spectral)
  .addBands(texture)
  .addBands(s2mean)
  .addBands([VV_mean, VH_mean])
  .clip(roi)
  .toFloat();

print('Final feature image:', allFeatures);


// ==========================================================
// 13. Sample extraction
// ==========================================================
var sample_features = allFeatures.sampleRegions({
  collection: samples,
  properties: ['class'],
  scale: 10,
  geometries: true
});

print('Sample feature table:', sample_features.limit(5));


// ==========================================================
// 14. Export
// ==========================================================
Export.table.toDrive({
  collection: sample_features,
  description: 'Sample_Features_' + year,
  folder: 'GEE_Export',
  fileNamePrefix: 'Sample_Features_' + year,
  fileFormat: 'CSV'
});
