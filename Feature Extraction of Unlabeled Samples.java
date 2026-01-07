/*************************************
 * Sentinel-2 NDVI Phenology & Feature Extraction (Annual)
 * Extract seasonal features, segmented integrals, spectral & texture features
 *************************************/

// --------------------- 1. Parameters ---------------------
var sp1 = table2;  // Sample points with class labels
var year = '2024';
var startDate = ee.Date(year + '-01-01');
var endDate = ee.Date(year + '-12-31');
var roi = table.geometry().bounds(); // Study area

// --------------------- 2. Load Sentinel-2 SR and calculate NDVI ---------------------
var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(function(img) {
    var ndvi = img.normalizedDifference(['B8','B4']).rename('NDVI');
    return ndvi.copyProperties(img, ['system:time_start']);
  });

// --------------------- 3. Monthly NDVI compositing ---------------------
var months = ee.List.sequence(1, 12);
var ndviMonthly = months.map(function(m) {
  var start = ee.Date(year + '-01-01').advance(ee.Number(m).subtract(1), 'month');
  var end = start.advance(1, 'month');
  var monthlyCol = s2.filterDate(start, end);
  var monthlyNDVI = monthlyCol.mean()
                    .set('month', m)
                    .set('system:time_start', start.millis());
  return monthlyNDVI;
});
var ndviCol = ee.ImageCollection(ndviMonthly);
var ndviList = ndviCol.toList(ndviCol.size());

// --------------------- 4. Determine growing season ---------------------
var ndviMax = ndviCol.max().rename('NDVImax');
var ndviThr = ndviMax.multiply(0.5).rename('NDVIthr');

var monthMaskImgs = months.map(function(m){
  m = ee.Number(m);
  var img = ee.Image(ndviList.get(m.subtract(1)));
  var mask = img.gt(ndviThr);
  return ee.Image.constant(m).updateMask(mask).toFloat();
});
var monthMaskCol = ee.ImageCollection.fromImages(monthMaskImgs);

var Month_start_val = monthMaskCol.min().unmask(0).rename('Month_start_val'); // SOS month
var Month_end_val = monthMaskCol.max().unmask(0).rename('Month_end_val');     // EOS month

// --------------------- 5. NDVI images within growing season ---------------------
var growingSeasonImgs = months.map(function(m){
  m = ee.Number(m);
  var nd = ee.Image(ndviList.get(m.subtract(1)));
  var cond = ee.Image.constant(m).gte(Month_start_val).and(ee.Image.constant(m).lte(Month_end_val));
  return nd.updateMask(cond).set('month', m);
});
var growingSeasonCol = ee.ImageCollection.fromImages(growingSeasonImgs);

// --------------------- 6. Phenology features ---------------------
// Start/End NDVI values
var startValImgs = months.map(function(m){
  m = ee.Number(m);
  var nd = ee.Image(ndviList.get(m.subtract(1)));
  return nd.updateMask(Month_start_val.eq(m));
});
var Start_val = ee.ImageCollection.fromImages(startValImgs).mosaic().unmask(0).rename('Start_val');

var endValImgs = months.map(function(m){
  m = ee.Number(m);
  var nd = ee.Image(ndviList.get(m.subtract(1)));
  return nd.updateMask(Month_end_val.eq(m));
});
var End_val = ee.ImageCollection.fromImages(endValImgs).mosaic().unmask(0).rename('End_val');

// Base value
var Base_val = Start_val.add(End_val).divide(2).rename('Base_val');

// Peak NDVI and amplitude
var Peak_val = growingSeasonCol.max().rename('Peak_val');
var Ampl = Peak_val.subtract(Base_val).rename('Ampl');

// --------------------- 7. Growth rate (derivative) ---------------------
var gsList = ee.List(months.map(function(m){
  m = ee.Number(m);
  return ee.Image(ndviList.get(m.subtract(1))).updateMask(
    ee.Image.constant(m).gte(Month_start_val).and(ee.Image.constant(m).lte(Month_end_val))
  );
}));

var diffImgs = ee.List.sequence(1,11).map(function(i){
  i = ee.Number(i);
  var before = ee.Image(gsList.get(i.subtract(1)));
  var after = ee.Image(gsList.get(i));
  return after.subtract(before).rename('diff').set('index', i);
});
var diffCol = ee.ImageCollection.fromImages(diffImgs);

var Max_increase = diffCol.max().rename('Max_increase');
var Min_decrease = diffCol.min().multiply(-1).rename('Min_decrease');

// --------------------- 8. Peak month ---------------------
var tol = 1e-6;
var peakMonthImgs = months.map(function(m){
  m = ee.Number(m);
  var nd = ee.Image(ndviList.get(m.subtract(1)));
  var cond = nd.subtract(Peak_val).abs().lte(tol);
  return ee.Image.constant(m).updateMask(cond).toFloat();
});
var peakMonthCol = ee.ImageCollection.fromImages(peakMonthImgs);
var Month_peak_val = peakMonthCol.min().unmask(0).rename('Month_peak_val');

// --------------------- 9. Segmented cumulative NDVI ---------------------
var prePeakImgs = months.map(function(m){
  m = ee.Number(m);
  var nd = ee.Image(ndviList.get(m.subtract(1)));
  var cond = ee.Image.constant(m).gte(Month_start_val).and(ee.Image.constant(m).lte(Month_peak_val));
  return nd.updateMask(cond);
});
var prePeakCol = ee.ImageCollection.fromImages(prePeakImgs);
var L_cumulative = prePeakCol.sum().unmask(0).rename('L_cumulative');

var postPeakImgs = months.map(function(m){
  m = ee.Number(m);
  var nd = ee.Image(ndviList.get(m.subtract(1)));
  var cond = ee.Image.constant(m).gte(Month_peak_val).and(ee.Image.constant(m).lte(Month_end_val));
  return nd.updateMask(cond);
});
var postPeakCol = ee.ImageCollection.fromImages(postPeakImgs);
var R_cumulative = postPeakCol.sum().unmask(0).rename('R_cumulative');

// --------------------- 10. Merge all phenology features ---------------------
var phenology = ndviMax.addBands([ndviThr, Month_start_val, Month_peak_val, Month_end_val,
                                  Start_val, End_val, Base_val, Peak_val, Ampl,
                                  L_cumulative, R_cumulative, Max_increase, Min_decrease
                                 ]).clip(roi).toFloat();

// --------------------- 11. Spectral indices (annual median) ---------------------
var s21 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(function(img){
    return img.divide(10000)
      .select(['B2','B3','B4','B8','B11'])
      .copyProperties(img,['system:time_start']);
  });

var s2mean = s21.median();
var ndvi = s2mean.normalizedDifference(['B8','B4']).rename('NDVI');
var ndwi = s2mean.normalizedDifference(['B3','B8']).rename('NDWI');
var gcvi = s2mean.expression('B8 / B3 - 1', {'B8': s2mean.select('B8'), 'B3': s2mean.select('B3')}).rename('GCVI');
var evi = s2mean.expression('2.5*(B8-B4)/(B8+6*B4-7.5*B2+1)', 
  {'B8': s2mean.select('B8'), 'B4': s2mean.select('B4'), 'B2': s2mean.select('B2')}).rename('EVI');
var savi = s2mean.expression('((B8-B4)/(B8+B4+0.5))*(1+0.5)', 
  {'B8': s2mean.select('B8'), 'B4': s2mean.select('B4')}).rename('SAVI');
var rvi = s2mean.expression('B8/B4', {'B8': s2mean.select('B8'), 'B4': s2mean.select('B4')}).rename('RVI');
var wi = s2mean.expression('(B3+2.5*B11)/(1.5*B8+2*B4+0.5*B11)',
  {'B3': s2mean.select('B3'),'B4': s2mean.select('B4'),'B8': s2mean.select('B8'),'B11': s2mean.select('B11')}).rename('WI');
var msi = s2mean.expression('B11/B8', {'B11': s2mean.select('B11'), 'B8': s2mean.select('B8')}).rename('MSI');
var ndbi = s2mean.normalizedDifference(['B11','B8']).rename('NDBI');
var ndmi = s2mean.normalizedDifference(['B8','B11']).rename('NDMI');
var mndwi = s2mean.normalizedDifference(['B3','B11']).rename('MNDWI');
var dvi = s2mean.expression('B8-B4', {'B8': s2mean.select('B8'),'B4': s2mean.select('B4')}).rename('DVI');

var spectral = ndvi.addBands([ndwi, gcvi, evi, savi, rvi, wi, msi, ndbi, ndmi, mndwi, dvi]);

// --------------------- 12. Texture features ---------------------
var gray = s2mean.expression('0.3*B8 + 0.59*B4 + 0.11*B3',
                             {'B8': s2mean.select('B8'),'B4': s2mean.select('B4'),'B3': s2mean.select('B3')}).rename('Gray');
var glcm = gray.multiply(100).toInt().glcmTexture({size: 3});
var texture = glcm.select(['Gray_contrast','Gray_asm','Gray_corr','Gray_var',
                           'Gray_idm','Gray_diss','Gray_sent','Gray_ent']);

// --------------------- 13. Sentinel-1 features ---------------------
var sentinel1 = ee.ImageCollection('COPERNICUS/S1_GRD')
                    .filterBounds(roi).filterDate(startDate,endDate)
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'))
                    .filter(ee.Filter.eq('instrumentMode','IW'));

var vvVhIwAsc = sentinel1.filter(ee.Filter.eq('orbitProperties_pass','ASCENDING'));
var vvVhIwDesc = sentinel1.filter(ee.Filter.eq('orbitProperties_pass','DESCENDING'));
var vvIwAscDescMean = vvVhIwAsc.merge(vvVhIwDesc).select('VV').mean().clip(roi);
var vhIwAscDescMean = vvVhIwAsc.merge(vvVhIwDesc).select('VH').mean().clip(roi);

// --------------------- 14. Merge all features ---------------------
var allFeatures = phenology.addBands(spectral)
                           .addBands(texture)
                           .addBands(s2mean)
                           .addBands(vvIwAscDescMean)
                           .addBands(vhIwAscDescMean)
                           .clip(roi)
                           .toFloat();

// --------------------- 15. Sampling & Export ---------------------
var image = allFeatures.select(['Peak_val','Base_val','B2','MNDWI','NDVI','RVI',
                                'End_val','B12','VH','B3','B4','VV','Start_val','B11','NDWI']);

var scale = 10;

// Sample the features at your sample points (keep original class property)
var sampledPoints = image.sampleRegions({
  collection: sp1,
  scale: scale,
  geometries: true
});
Export.table.toDrive({
  collection: sampledPoints,
  description: 'Samples_With_Features',
  folder: 'GEE_Export',
  fileNamePrefix: 'samples_with_features',
  fileFormat: 'CSV'
});

// Batch random sampling for large dataset
var batchSize = 5000;
for (var i = 1; i <= 6; i++) {
  var samples = image.sample({
    region: roi,
    scale: scale,
    numPixels: batchSize,
    geometries: true,
    seed: i*1234
  });
  Export.table.toDrive({
    collection: samples,
    description: 'Sentinel_Samples_seed_' + i,
    folder: 'GEE_Samples',
    fileNamePrefix: 'sentinel_samples_seed_' + i,
    fileFormat: 'CSV'
  });
}

