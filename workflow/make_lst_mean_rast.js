// Define the bounding box.
var bbox = ee.Geometry.Rectangle([10, 47, 15, 50]);

// Load the Daymet dataset.
var daymet = ee.ImageCollection("MODIS/061/MOD11A1").filterDate('2000-08-30', '2010-08-30');


// Get the average tmax for the last 30 years.
var tmaxMean = daymet.mean();

// Convert the ImageCollection to an Image.
//var tmaxMeanImage = ee.Image(tmaxMean);

// Clip the tmaxMean image to the bounding box.
//var clipped = tmaxMeanImage.clip(bbox);
var clipped = tmaxMean.clip(bbox);

// Define the visualization parameters.
var visParams = {
  bands: "LST_Day_1km",
  min: -10,
  max: 40,
  palette: ['blue', 'white', 'red']
};

// Add the tmaxMean image to the map.
Map.addLayer(clipped, visParams, 'Mean MODIS LST');

// Center and zoom to the bounding box.
Map.centerObject(bbox);
Map.setZoom(9);

// Define the export parameters.
var exportParams = {
  image: clipped,
  description: 'modis_lst_mean',
  region: bbox,
  scale: 1000, // set the export resolution
  crs: 'EPSG:4326' // set the export projection
};

// Export the image to Google Drive.
Export.image.toDrive(exportParams);