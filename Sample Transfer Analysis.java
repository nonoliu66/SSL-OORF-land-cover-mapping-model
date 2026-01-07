/*****************************************************
 * CLCD-based Multi-year Change Detection
 * + Convolution-based Pure Pixel Screening
 * + Unchanged-area Stratified Sampling
 * + Sensor-consistent Spectral Similarity Analysis
 *****************************************************/

// ----------------------------------------------------
// 1. Input data
// ----------------------------------------------------
var roi = table.geometry().bounds(); // Study area

var clcd2019 = ee.Image(image);
var clcd2020 = ee.Image(image2);
var clcd2021 = ee.Image(image3);
var clcd2022 = ee.Image(image4);
var clcd2023 = ee.Image(image5);
var clcd2024 = ee.Image(image6);

// ----------------------------------------------------
// 2. Preprocessing: band selection and clipping
// ----------------------------------------------------
function prep(img) {
  return img.select('b1').clip(roi).toInt();
}

clcd2019 = prep(clcd2019);
clcd2020 = prep(clcd2020);
clcd2021 = prep(clcd2021);
clcd2022 = prep(clcd2022);
clcd2023 = prep(clcd2023);
clcd2024 = prep(clcd2024);

// ----------------------------------------------------
// 3. Common spatial mask (valid pixels in all years)
// ----------------------------------------------------
var commonMask = clcd2019.addBands(clcd2020)
  .addBands(clcd2021)
  .addBands(clcd2022)
  .addBands(clcd2023)
  .addBands(clcd2024)
  .reduce(ee.Reducer.min())
  .gt(0)
  .rename('mask');

clcd2019 = clcd2019.updateMask(commonMask);
clcd2020 = clcd2020.updateMask(commonMask);
clcd2021 = clcd2021.updateMask(commonMask);
clcd2022 = clcd2022.updateMask(commonMask);
clcd2023 = clcd2023.updateMask(commonMask);
clcd2024 = clcd2024.updateMask(commonMask);

// ----------------------------------------------------
// 4. Multi-year change detection
// ----------------------------------------------------
var minImg = clcd2019.addBands(clcd2020)
  .addBands(clcd2021)
  .addBands(clcd2022)
  .addBands(clcd2023)
  .addBands(clcd2024)
  .reduce(ee.Reducer.min())
  .rename('min');

var maxImg = clcd2019.addBands(clcd2020)
  .addBands(clcd2021)
  .addBands(clcd2022)
  .addBands(clcd2023)
  .addBands(clcd2024)
  .reduce(ee.Reducer.max())
  .rename('max');

// Pixels unchanged across all years
var unchanged = minImg.eq(maxImg).selfMask().rename('unchanged');
// Pixels experiencing changes
var changed = minImg.neq(maxImg).selfMask().rename('changed');

// ----------------------------------------------------
// 5. Convolution-based pure pixel screening
// ----------------------------------------------------
var kernel = ee.Kernel.square(1); // 3×3 neighborhood

var unchangedConv = unchanged.unmask(0).toFloat().convolve(kernel);
var diffUnchanged = unchanged.toFloat().subtract(unchangedConv).abs();

// Retain only spatially homogeneous unchanged pixels
var pureUnchanged = unchanged.updateMask(diffUnchanged.eq(0));

// ----------------------------------------------------
// 6. Sample extraction from unchanged & pure areas
// ----------------------------------------------------
var pureClass = clcd2023
  .updateMask(pureUnchanged)
  .rename('b1')
  .clip(roi);

var samples = pureClass.stratifiedSample({
  numPoints: 1000,
  classBand: 'b1',
  region: roi,
  scale: 30,
  geometries: true,
  seed: 2025
});

Map.addLayer(samples, {color: 'red'}, 'Pure unchanged samples');

// ----------------------------------------------------
// 7. Landsat preprocessing functions
// ----------------------------------------------------
var rmCloud = function(image) {
  var qaMask = image.select('QA_PIXEL')
    .bitwiseAnd(parseInt('11111', 2))
    .eq(0);
  var satMask = image.select('QA_RADSAT').eq(0);
  return image.updateMask(qaMask).updateMask(satMask);
};

function applyScaleFactors(image) {
  var optical = image.select('SR_B.')
    .multiply(0.0000275)
    .add(-0.2);
  var thermal = image.select('ST_B.*')
    .multiply(0.00341802)
    .add(149.0);
  return image.addBands(optical, null, true)
              .addBands(thermal, null, true);
}

// Band definitions
var L8 = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'];
var L7 = ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7'];
var Bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2'];

// ----------------------------------------------------
// 8. Reference year imagery (2023): Landsat 7 + 8
// ----------------------------------------------------
var Year = 2023;

var LS8_2023 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filter(ee.Filter.calendarRange(Year, Year, 'year'))
  .map(rmCloud)
  .map(applyScaleFactors)
  .select(L8, Bands);

var LS7_2023 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterBounds(roi)
  .filter(ee.Filter.calendarRange(Year, Year, 'year'))
  .map(rmCloud)
  .map(applyScaleFactors)
  .select(L7, Bands);

var refImage = LS8_2023.merge(LS7_2023)
  .mean()
  .clip(roi);

// ----------------------------------------------------
// 9. Target year imagery (2019): Landsat 7 + 8
// ----------------------------------------------------
Year = 2019;

var LS8_2019 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filter(ee.Filter.calendarRange(Year, Year, 'year'))
  .map(rmCloud)
  .map(applyScaleFactors)
  .select(L8, Bands);

var LS7_2019 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterBounds(roi)
  .filter(ee.Filter.calendarRange(Year, Year, 'year'))
  .map(rmCloud)
  .map(applyScaleFactors)
  .select(L7, Bands);

var targetImage = LS8_2019.merge(LS7_2019)
  .mean()
  .clip(roi);

// ----------------------------------------------------
// 10. Spectral similarity computation (SAD & ED)
// ----------------------------------------------------
var sadImg = refImage
  .spectralDistance(targetImage, 'sam')
  .cos()
  .rename('SAD');

var edImg = refImage
  .spectralDistance(targetImage, 'sed')
  .rename('ED');

// ----------------------------------------------------
// 11. Extract spectral similarity values at sample points
// ----------------------------------------------------
var samplesFinal = sadImg
  .addBands(edImg)
  .addBands(targetImage)
  .sampleRegions({
    collection: samples,
    properties: ['b1'],
    scale: 30,
    geometries: true
  })
  .map(function(f) {
    var c = f.geometry().coordinates();
    return f.set({
      longitude: c.get(0),
      latitude: c.get(1)
    });
  });

// ----------------------------------------------------
// 12. Export results
// ----------------------------------------------------
Export.table.toDrive({
  collection: samplesFinal,
  description: 'SamplePoints_2019_2023_L7_L8',
  fileFormat: 'CSV'
});
