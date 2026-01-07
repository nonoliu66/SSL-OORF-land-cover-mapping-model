// ======================= Parameters =======================
var startDate = '2023-01-01';
var endDate   = '2023-12-31';
var numClusters = 30;                      // Number of clusters
var region = table.geometry().bounds();    // Study area

// =================== Sentinel-2 Preprocessing ===================
// Cloud and shadow masking using SCL band
// Retained classes: vegetation (4), bare soil (5), water (6), snow/ice (11)
function maskS2(image) {
  var scl = image.select('SCL');
  var mask = scl.eq(4)
                .or(scl.eq(5))
                .or(scl.eq(6))
                .or(scl.eq(11));

  return image.updateMask(mask)
    .select([
      'B2',   // Blue
      'B3',   // Green
      'B4',   // Red
      'B8',   // NIR
      'B8A',  // Narrow NIR
      'B11',  // SWIR-1
      'B12'   // SWIR-2
    ])
    .divide(10000)                         // Convert DN to surface reflectance
    .copyProperties(image, ['system:time_start']);
}

// Load Sentinel-2 L2A imagery
var s2Col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate(startDate, endDate)
  .filterBounds(region)
  .map(maskS2);

// Annual median composite
var s2Composite = s2Col.median().clip(region);

// =================== Sentinel-1 Preprocessing ===================
// Load Sentinel-1 GRD imagery (IW mode, VV and VH polarizations)
var s1Col = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterDate(startDate, endDate)
  .filterBounds(region)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('resolution_meters', 10));

// Annual median composite (linear backscatter)
var s1Median = s1Col.select(['VV', 'VH']).median();

// Convert backscatter to decibel (dB)
var s1dB = s1Median.log10().multiply(10)
  .rename(['VV_dB', 'VH_dB']);

// =================== Feature Fusion ===================
var fusedImage = s2Composite.addBands(s1dB);

// =================== Unsupervised Clustering ===================
// Randomly sample pixels from the fused feature space
var trainingSamples = fusedImage.sample({
  region: region,
  scale: 10,              // Match Sentinel spatial resolution
  numPixels: 5000,
  geometries: false
});

// Train KMeans clusterer
var clusterer = ee.Clusterer.wekaKMeans(numClusters)
  .train(trainingSamples);

// Apply clustering to the fused image
var clusteredImage = fusedImage.cluster(clusterer);

// =================== Visualization ===================
Map.addLayer(
  s2Composite,
  {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3},
  'Sentinel-2 RGB Composite'
);

Map.addLayer(
  clusteredImage.randomVisualizer(),
  {},
  'KMeans Clustering Result'
);

// =================== Export ===================
Export.image.toDrive({
  image: clusteredImage,
  description: 'unsupervised_classification_30classes_2023',
  folder: 'EarthEngineExports',
  fileNamePrefix: 'unsupervised_classification_30classes_2023',
  region: region,
  scale: 10,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
