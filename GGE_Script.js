// The Google Earth Engine (GEE) source code and data used in:
// "Two Decades Mangroves Loss Monitoring Using Random Forest and Landsat Data in East Luwu, Indonesia (2000–2020)"
// published in Geomatics MDPI (2022) https://doi.org/10.3390/geomatics2030016

//Copy from: https://code.earthengine.google.com/8e5aff637aee53ecf84e5c0aefa90623

//Inspired by: Barenblitt, A.; Fatoyinbo, T. (2020). Remote Sensing for Mangroves in Support of the UN Sustainable Development Goals. NASA Applied Remote Sensing Training Program (ARSET).
//https://appliedsciences.nasa.gov/get-involved/training/english/arset-remote-sensing-mangroves-support-un-sustainable-development

var East_Luwu = ee.FeatureCollection("users/SINAUGIS/East_Luwu")
Map.centerObject(coastal,10)
Map.addLayer(East_Luwu)

//Cloud masking funtion for Landsat-8
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask).divide(10000);
}

//Cloud masking funtion for Landsat-8
var cloudMaskL457 = function(image) {
  var qa = image.select('pixel_qa');
  // If the cloud bit (5) is set and the cloud confidence (7) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var cloud = qa.bitwiseAnd(1 << 5)
                  .and(qa.bitwiseAnd(1 << 7))
                  .or(qa.bitwiseAnd(1 << 3));
  // Remove edge pixels that don't occur in all bands
  var mask2 = image.mask().reduce(ee.Reducer.min());
  return image.updateMask(cloud.not()).updateMask(mask2).divide(10000);
};

////////////-----------------------------------------------------//////////
//////////// -- Landsat-7 Processing for the year 2000 - 2013 -- //////////
////////////-----------------------------------------------------//////////
var addIndicesL7 = function(img) {
  // NDVI
  var ndvi = img.normalizedDifference(['B4','B3']).rename('NDVI');
  
  // MNDWI (Modified Normalized Difference Water Index - Hanqiu Xu, 2006)
  var mndwi = img.normalizedDifference(['B2','B5']).rename('MNDWI');
  
  //EVI
  var cmri = img.expression('NDVI-NDWI',{
    'NDVI':img.normalizedDifference(['B4','B3']),
    'NDWI':img.normalizedDifference(['B2','B4'])
  }).rename('CMRI');
  
  // NDMI (Normalized Difference Mangrove Index - Shi et al 2016 - New spectral metrics for mangrove forest identification)
  var ndmi = img.normalizedDifference(['B7','B2']).rename('NDMI');
   
  // MMRI
  var mmri = img.expression('abs((GREEN-MIR)/(GREEN+MIR))-abs((NIR-RED)/(NIR+RED))/abs((GREEN-MIR)/(GREEN+MIR))+abs((NIR-RED)/(NIR+RED))',{
    'NIR':img.select('B4'),
    'GREEN':img.select('B2'),
    'MIR':img.select('B5'),
    'RED':img.select('B3')
  }).rename('MMRI');
  
  var blue = img.select('B1').rename('BLUE')
  var green = img.select('B2').rename('GREEN')
  var red = img.select('B3').rename('RED')
  var nir = img.select('B4').rename('NIR')
  var swir1 = img.select('B5').rename('SWIR1')
  var swir2 = img.select('B7').rename('SWIR2')
  
  
  return img
    .addBands(blue)
    .addBands(green)
    .addBands(red)
    .addBands(nir)
    .addBands(swir1)
    .addBands(swir2)
    .addBands(ndvi)
    .addBands(cmri)
    .addBands(ndmi)
    .addBands(mmri)
    .addBands(mndwi);
};

var year00 = 2000; 
var startDate00 = (year00)+'-01-01'; 
var endDate00 = (year00)+'-12-31'; 

var l700 = L7.filterDate(startDate00,endDate00)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite00 = l700
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var srtmClip = SRTM.clip(East_Luwu);
var elevationMask = srtmClip.lt(40);
var MNDWIMask00 = composite00.select('MNDWI').lt(0);
var compositeNew00 = composite00
                        .updateMask(MNDWIMask00)
                        .updateMask(elevationMask)
                        
var visPar = {bands:['B4','B5','B3'], min: 0, max: 0.35}; 
print(compositeNew00)
Map.addLayer(compositeNew00.clip(East_Luwu), visPar, 'Landsat Composite')

///Build Sample training sample for L7 by using the year 2000
var classes = Mangrove.merge(NonMangrove)
Export.table.toDrive({collection: classes,
                      description: 'Training_Landsat7',
                      fileFormat: 'SHP'});
                      
var bands = ['BLUE','GREEN','RED','NIR','SWIR1','SWIR2','NDVI','CMRI','NDMI','MMRI']
var image00 = compositeNew00.select(bands).clip(East_Luwu)
   
var samples = image00.sampleRegions({
    collection: classes,
    properties: ['landcover'],
    scale: 30
    }).randomColumn('random');

var TotMangrove = image00.sampleRegions({
    collection: Mangrove,
    properties: ['landcover'],
    scale: 30
    }).randomColumn('random');

var TotNonMangrove = image00.sampleRegions({
    collection: NonMangrove,
    properties: ['landcover'],
    scale: 30
    }).randomColumn('random');


print('Samples n =', samples.aggregate_count('.all'));
print('Samples Mangrove =', TotMangrove.aggregate_count('.all'))
print('Samples NonMangrove =', TotNonMangrove.aggregate_count('.all'))


///Randfom forest model for Landat-7 
var classifier = ee.Classifier.smileRandomForest(500).train({ 
    features: samples.select(['BLUE','GREEN','RED','NIR','SWIR1','SWIR2','NDVI','CMRI','NDMI','MMRI', 'landcover']),
    classProperty: 'landcover',
    inputProperties: bands
    });
print('Classifier_L7:',classifier.explain()); 


///Applied Classifier for the year 2000
var classifiedrf00 = image00.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount00 = classifiedrf00.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask00 = pixelcount00.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask00 = classifiedrf00.select('classification').gt(0)
var classed00= classifiedrf00.updateMask(countmask00).updateMask(classMask00).clip(coastal)

Map.addLayer (classed00, {min: 1, max: 1, palette:'blue'}, 'Mangrove Extent 2000');
var imageex00 = composite00.select('B4','B5','B3');

var mangrove = ee.ImageCollection('LANDSAT/MANGROVE_FORESTS');
var mangrovesVis = {
  min: 0,
  max: 1.0,
  palette: ['d40115'],
};
Map.addLayer(mangrove.first().clip(coastal), mangrovesVis, 'Mangroves Giri');

////////////////////////////////////////
///////////Applied 2001//////////////////
////////////////////////////////////////

var year01 = 2001;
var startDate01 = (year01)+'-01-01'; 
var endDate01 = (year01)+'-12-31'; 

var l701 = L7.filterDate(startDate01,endDate01)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite01 = l701
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask01 = composite01.select('MNDWI').lt(0);
var compositeNew01 = composite01
                        .updateMask(MNDWIMask01)
                        .updateMask(elevationMask);

// Map.addLayer(compositeNew01.clip(East_Luwu), visPar, 'Landsat Composite 2001')


///APLLIED MODEL to 2001
var image01 = compositeNew01.select(bands).clip(East_Luwu)
var classifiedrf01 = image01.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount01 = classifiedrf01.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask01 = pixelcount01.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask01 = classifiedrf01.select('classification').gt(0)
var classed01= classifiedrf01.updateMask(countmask01).updateMask(classMask01).clip(coastal)
// Map.addLayer (classed01, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2001');
var imageex01 = composite01.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2002//////////////////
////////////////////////////////////////

var year02 = 2002;
var startDate02 = (year02)+'-01-01'; 
var endDate02 = (year02)+'-12-31'; 

var l702 = L7.filterDate(startDate02,endDate02)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite02 = l702
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask02 = composite02.select('MNDWI').lt(0);
var compositeNew02 = composite02
                        .updateMask(MNDWIMask02)
                        .updateMask(elevationMask);

// Map.addLayer(compositeNew02.clip(East_Luwu), visPar, 'Landsat Composite 2002')


///APLLIED MODEL to 2002
var image02 = compositeNew02.select(bands).clip(East_Luwu)
    var classifiedrf02 = image02.select(bands)
                      .classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount02 = classifiedrf02.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask02 = pixelcount02.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask02 = classifiedrf02.select('classification').gt(0)
var classed02= classifiedrf02.updateMask(countmask02).updateMask(classMask02).clip(coastal)
// Map.addLayer (classed02, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2002');
var imageex02 = composite02.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2003//////////////////
////////////////////////////////////////

var year03 = 2003;
var startDate03 = (year03)+'-01-01'; 
var endDate03 = (year03)+'-12-31'; 

var l703 = L7.filterDate(startDate03,endDate03)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite03 = l703
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask03 = composite03.select('MNDWI').lt(0);
var compositeNew03 = composite03
                        .updateMask(MNDWIMask03)
                        .updateMask(elevationMask);

// Map.addLayer(compositeNew03.clip(East_Luwu), visPar, 'Landsat Composite 2003')


///APLLIED MODEL to 2003
var image03 = compositeNew03.select(bands).clip(East_Luwu)
var classifiedrf03 = image03.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount03 = classifiedrf03.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask03 = pixelcount03.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask03 = classifiedrf03.select('classification').gt(0)
var classed03= classifiedrf03.updateMask(countmask03).updateMask(classMask03).clip(coastal)
// Map.addLayer (classed03, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2003');
var imageex03 = composite03.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2004//////////////////
////////////////////////////////////////

var year04 = 2004;
var startDate04 = (year04)+'-01-01'; 
var endDate04 = (year04)+'-12-31'; 

var l704 = L7.filterDate(startDate04,endDate04)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite04 = l704
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask04 = composite04.select('MNDWI').lt(0);
var compositeNew04 = composite04
                        .updateMask(MNDWIMask04)
                        .updateMask(elevationMask);

// Map.addLayer(compositeNew04.clip(East_Luwu), visPar, 'Landsat Composite 2004')


///APLLIED MODEL to 2004
var image04 = compositeNew04.select(bands).clip(East_Luwu)
var classifiedrf04 = image04.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount04 = classifiedrf04.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask04 = pixelcount04.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask04 = classifiedrf04.select('classification').gt(0)
var classed04= classifiedrf04.updateMask(countmask04).updateMask(classMask04).clip(coastal)
// Map.addLayer (classed04, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2004');
var imageex04 = composite04.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2005//////////////////
////////////////////////////////////////

var year05 = 2005;
var startDate05 = (year05)+'-01-01'; 
var endDate05 = (year05)+'-12-31'; 

var l705 = L7.filterDate(startDate05,endDate05)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite05 = l705
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask05 = composite05.select('MNDWI').lt(0);
var compositeNew05 = composite05
                        .updateMask(MNDWIMask05)
                        .updateMask(elevationMask);

Map.addLayer(compositeNew05.clip(East_Luwu), visPar, 'Landsat Composite 2005')


///APLLIED MODEL to 2005
var image05 = compositeNew05.select(bands).clip(East_Luwu)
var classifiedrf05 = image05.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount05 = classifiedrf05.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask05 = pixelcount05.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask05 = classifiedrf05.select('classification').gt(0)
var classed05= classifiedrf05.updateMask(countmask05).updateMask(classMask05).clip(coastal)
Map.addLayer (classed05, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2005');
var imageex05 = composite05.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2006//////////////////
////////////////////////////////////////

var year06 = 2006;
var startDate06 = (year06)+'-01-01'; 
var endDate06 = (year06)+'-12-31'; 

var l706 = L7.filterDate(startDate06,endDate06)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite06 = l706
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask06 = composite06.select('MNDWI').lt(0);
var compositeNew06 = composite06
                        .updateMask(MNDWIMask06)
                        .updateMask(elevationMask);

// Map.addLayer(compositeNew06.clip(East_Luwu), visPar, 'Landsat Composite 2006')


///APLLIED MODEL to 2006
var image06 = compositeNew06.select(bands).clip(East_Luwu)
var classifiedrf06 = image06.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount06 = classifiedrf06.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask06 = pixelcount06.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask06 = classifiedrf06.select('classification').gt(0)
var classed06= classifiedrf06.updateMask(countmask06).updateMask(classMask06).clip(coastal)
// Map.addLayer (classed06, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2006');
var imageex06 = composite06.select('B4','B5','B3');


////////////////////////////////////////
///////////Applied 2007//////////////////
////////////////////////////////////////

var year07 = 2007;
var startDate07 = (year07)+'-01-01'; 
var endDate07 = (year07)+'-12-31'; 

var l707 = L7.filterDate(startDate07,endDate07)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite07 = l707
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask07 = composite07.select('MNDWI').lt(0);
var compositeNew07 = composite07
                        .updateMask(MNDWIMask07)
                        .updateMask(elevationMask);

// Map.addLayer(compositeNew07.clip(East_Luwu), visPar, 'Landsat Composite 2007')


///APLLIED MODEL to 2007
var image07 = compositeNew07.select(bands).clip(East_Luwu)
var classifiedrf07 = image07.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount07 = classifiedrf07.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask07 = pixelcount07.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask07 = classifiedrf07.select('classification').gt(0)
var classed07= classifiedrf07.updateMask(countmask07).updateMask(classMask07).clip(coastal)
// Map.addLayer (classed07, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2007');
var imageex07 = composite07.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2008//////////////////
////////////////////////////////////////

var year08 = 2008;
var startDate08 = (year08)+'-01-01'; 
var endDate08 = (year08)+'-12-31'; 

var l708 = L7.filterDate(startDate08,endDate08)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite08 = l708
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask08 = composite08.select('MNDWI').lt(0);
var compositeNew08 = composite08
                        .updateMask(MNDWIMask08)
                        .updateMask(elevationMask);

// Map.addLayer(compositeNew07.clip(East_Luwu), visPar, 'Landsat Composite 2008')


///APLLIED MODEL to 2008
var image08 = compositeNew08.select(bands).clip(East_Luwu)
var classifiedrf08 = image08.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount08 = classifiedrf08.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask08 = pixelcount08.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask08 = classifiedrf08.select('classification').gt(0)
var classed08= classifiedrf08.updateMask(countmask08).updateMask(classMask08).clip(coastal)
// Map.addLayer (classed08, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2008');
var imageex08 = composite08.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2009//////////////////
////////////////////////////////////////

var year09 = 2009;
var startDate09 = (year09)+'-01-01'; 
var endDate09 = (year09)+'-12-31'; 

var l709 = L7.filterDate(startDate09,endDate09)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite09 = l709
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask09 = composite09.select('MNDWI').lt(0);
var compositeNew09 = composite09
                        .updateMask(MNDWIMask09)
                        .updateMask(elevationMask);
// Map.addLayer(compositeNew07.clip(East_Luwu), visPar, 'Landsat Composite 2009')


///APLLIED MODEL to 2009
var image09 = compositeNew09.select(bands).clip(East_Luwu)
var classifiedrf09 = image09.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount09 = classifiedrf09.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask09 = pixelcount09.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask09 = classifiedrf09.select('classification').gt(0)
var classed09= classifiedrf09.updateMask(countmask09).updateMask(classMask09).clip(coastal)
// Map.addLayer (classed09, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2009');
var imageex09 = composite09.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2010//////////////////
////////////////////////////////////////

var year10 = 2010;
var startDate10 = (year10)+'-01-01'; 
var endDate10 = (year10)+'-12-31'; 

var l710 = L7.filterDate(startDate10,endDate10)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite10 = l710
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask10 = composite10.select('MNDWI').lt(0);
var compositeNew10 = composite10
                        .updateMask(MNDWIMask10)
                        .updateMask(elevationMask);

Map.addLayer(compositeNew10.clip(East_Luwu), visPar, 'Landsat Composite 2010')


///APLLIED MODEL to 2010
var image10 = compositeNew10.select(bands).clip(East_Luwu)
var classifiedrf10 = image10.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount10 = classifiedrf10.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask10 = pixelcount10.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask10 = classifiedrf10.select('classification').gt(0)
var classed10= classifiedrf10.updateMask(countmask10).updateMask(classMask10).clip(coastal)
Map.addLayer (classed10, {min: 1, max: 1, palette:'purple'}, 'Mangrove Extent 2010');
var imageex10 = composite10.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2011//////////////////
////////////////////////////////////////

var year11 = 2011;
var startDate11 = (year11)+'-01-01'; 
var endDate11 = (year11)+'-12-31'; 

var l711 = L7.filterDate(startDate11,endDate11)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite11 = l711
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask11 = composite11.select('MNDWI').lt(0);
var compositeNew11 = composite11
                        .updateMask(MNDWIMask11)
                        .updateMask(elevationMask);

// Map.addLayer(compositeNew11.clip(East_Luwu), visPar, 'Landsat Composite 2011')


///APLLIED MODEL to 2011
var image11 = compositeNew11.select(bands).clip(East_Luwu)
var classifiedrf11 = image11.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount11 = classifiedrf11.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask11 = pixelcount11.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask11 = classifiedrf11.select('classification').gt(0)
var classed11= classifiedrf11.updateMask(countmask11).updateMask(classMask11).clip(coastal)
// Map.addLayer (classed11, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2011');
var imageex11 = composite11.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2012//////////////////
////////////////////////////////////////

var year12 = 2012;
var startDate12 = (year12)+'-01-01'; 
var endDate12 = (year12)+'-12-31'; 

var l712 = L7.filterDate(startDate12,endDate12)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite12 = l712
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask12 = composite12.select('MNDWI').lt(0);
var compositeNew12 = composite12
                        .updateMask(MNDWIMask12)
                        .updateMask(elevationMask);

// Map.addLayer(compositeNew12.clip(East_Luwu), visPar, 'Landsat Composite 2012')


///APLLIED MODEL to 2012
var image12 = compositeNew12.select(bands).clip(East_Luwu)
var classifiedrf12 = image12.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount12 = classifiedrf12.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask12 = pixelcount12.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask12 = classifiedrf12.select('classification').gt(0)
var classed12= classifiedrf12.updateMask(countmask12).updateMask(classMask12).clip(coastal)
// Map.addLayer (classed12, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2012');
var imageex12 = composite12.select('B4','B5','B3');

////////////////////////////////////////
///////////Applied 2013//////////////////
////////////////////////////////////////

var year13 = 2013;
var startDate13 = (year13)+'-01-01'; 
var endDate13 = (year13)+'-12-31'; 

var l713 = L7.filterDate(startDate13,endDate13)
    .map(cloudMaskL457)
    .map(addIndicesL7)

var composite13 = l713
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask13 = composite13.select('MNDWI').lt(0);
var compositeNew13 = composite13
                        .updateMask(MNDWIMask13)
                        .updateMask(elevationMask);

// Map.addLayer(compositeNew12.clip(East_Luwu), visPar, 'Landsat Composite 2013')/////////////////////CHANGE


///APLLIED MODEL to 2013
var image13 = compositeNew13.select(bands).clip(East_Luwu)
var classifiedrf13 = image13.select(bands).classify(classifier); 
                      
///pixel count to reduc noise
var pixelcount13 = classifiedrf13.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask13 = pixelcount13.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask13 = classifiedrf13.select('classification').gt(0)
var classed13= classifiedrf13.updateMask(countmask13).updateMask(classMask13).clip(coastal)
// Map.addLayer (classed13, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2013');
var imageex13 = composite13.select('B4','B5','B3');

/////////////////////////////// -- End of Landsat-7 Processing -- ///////////////////////////////


////////////-----------------------------------------------------//////////
//////////// -- Landsat-8 Processing for the year 2014 - 2020 -- //////////
////////////-----------------------------------------------------//////////
var addIndicesL8 = function(img) {
  // NDVI
  var ndvi = img.normalizedDifference(['B5','B4']).rename('NDVI');
  
  // MNDWI (Modified Normalized Difference Water Index - Hanqiu Xu, 2006)
  var mndwi = img.normalizedDifference(['B3','B6']).rename('MNDWI');
  
  //CMRI
  var cmri = img.expression('NDVI-NDWI',{
    'NDVI':img.normalizedDifference(['B5','B4']),
    'NDWI':img.normalizedDifference(['B3','B5']),
  }).rename('CMRI');
  
  // NDMI (Normalized Difference Mangrove Index - Shi et al 2016 - New spectral metrics for mangrove forest identification)
  var ndmi = img.normalizedDifference(['B7','B3']).rename('NDMI');
   
  // MMRI
  var mmri = img.expression('abs((GREEN-MIR)/(GREEN+MIR))-abs((NIR-RED)/(NIR+RED))/abs((GREEN-MIR)/(GREEN+MIR))+abs((NIR-RED)/(NIR+RED))',{
    'NIR':img.select('B5'),
    'GREEN':img.select('B3'),
    'MIR':img.select('B6'),
    'RED':img.select('B4')
  }).rename('MMRI');
  
  var blue = img.select('B2').rename('BLUE')
  var green = img.select('B3').rename('GREEN')
  var red = img.select('B4').rename('RED')
  var nir = img.select('B5').rename('NIR')
  var swir1 = img.select('B6').rename('SWIR1')
  var swir2 = img.select('B7').rename('SWIR2')
  
  
  return img
    .addBands(blue)
    .addBands(green)
    .addBands(red)
    .addBands(nir)
    .addBands(swir1)
    .addBands(swir2)
    .addBands(ndvi)
    .addBands(cmri)
    .addBands(ndmi)
    .addBands(mmri)
    .addBands(mndwi);
};

var visParL8 = {bands:['B5','B6','B4'], min: 0, max: 0.35};

////////////////////////////////////////
///////////Applied 2014//////////////////
////////////////////////////////////////

var year14 = 2014;
var startDate14 = (year14)+'-01-01'; 
var endDate14 = (year14)+'-12-31'; 

var l814 = L8.filterDate(startDate14,endDate14)
    .map(maskL8sr)
    .map(addIndicesL8)

var composite14 = l814
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask14 = composite14.select('MNDWI').lt(0.07);
var compositeNew14 = composite14
                        .updateMask(MNDWIMask14)
                        .updateMask(elevationMask);
                        
Map.addLayer(compositeNew14.clip(East_Luwu), visParL8, 'Landsat 8 Composite 2014')

///Build Sample training sample for L8 by using the year 2014
var classes_8 = mangrove_8.merge(nonmangrove_8)
Export.table.toDrive({collection: classes_8,
                      description: 'Training_Landsat8',
                      fileFormat: 'SHP'});

var bands = ['BLUE','GREEN','RED','NIR','SWIR1','SWIR2','NDVI','CMRI','NDMI','MMRI']
var image14 = compositeNew14.select(bands).clip(East_Luwu)
   
var samples_8 = image14.sampleRegions({
    collection: classes_8,
    properties: ['landcover'],
    scale: 30
    }).randomColumn('random');

var TotMangrove_8  = image14.sampleRegions({
    collection: mangrove_8,
    properties: ['landcover'],
    scale: 30
    }).randomColumn('random');

var TotNonMangrove_8  = image14.sampleRegions({
    collection: nonmangrove_8,
    properties: ['landcover'],
    scale: 30
    }).randomColumn('random');

    print('Samples_8  n =', samples_8.aggregate_count('.all'));
    print('Samples_8  Mangrove =', TotMangrove_8.aggregate_count('.all'))
    print('Samples_8  NonMangrove =', TotNonMangrove_8.aggregate_count('.all'))

///Randfom Forest model Landsat-8
    var classifier_8  = ee.Classifier.smileRandomForest(500).train({ 
    features: samples_8 .select(['BLUE','GREEN','RED','NIR','SWIR1','SWIR2','NDVI','CMRI','NDMI','MMRI', 'landcover']),
    classProperty: 'landcover',
    inputProperties: bands
    });

print('Classifier_L8:',classifier_8.explain()); 

///APLLIED MODEL to 2014
var classifiedrf14 = image14.select(bands).classify(classifier_8); 
                      
///pixel count to reduc noise
var pixelcount14 = classifiedrf14.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask14 = pixelcount14.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask14 = classifiedrf14.select('classification').gt(0)
var classed14= classifiedrf14.updateMask(countmask14).updateMask(classMask14).clip(coastal)
Map.addLayer (classed14, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2014');
var imageex14 = composite14.select('B5','B6','B4');

////////////////////////////////////////
///////////Applied 2015//////////////////
////////////////////////////////////////

var year15 = 2015;
var startDate15 = (year15)+'-01-01'; 
var endDate15 = (year15)+'-12-31'; 

var l815 = L8.filterDate(startDate15,endDate15)
    .map(maskL8sr)
    .map(addIndicesL8)

var composite15 = l815
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask15 = composite15.select('MNDWI').lt(0.07);
var compositeNew15 = composite15
                        .updateMask(MNDWIMask15)
                        .updateMask(elevationMask);
                        
Map.addLayer(compositeNew15.clip(East_Luwu), visParL8, 'Landsat 8 Composite 2015')


///APLLIED MODEL to 2015
var image15 = compositeNew15.select(bands).clip(East_Luwu)
var classifiedrf15 = image15.select(bands).classify(classifier_8); 
                      
///pixel count to reduc noise
var pixelcount15 = classifiedrf15.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask15 = pixelcount15.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask15 = classifiedrf15.select('classification').gt(0)
var classed15= classifiedrf15.updateMask(countmask15).updateMask(classMask15)
Map.addLayer (classed15, {min: 1, max: 1, palette:'brown'}, 'Mangrove Extent 2015');
var imageex15 = composite15.select('B5','B6','B4');

////////////////////////////////////////
///////////Applied 2016//////////////////
////////////////////////////////////////

var year16 = 2016;
var startDate16 = (year16)+'-01-01'; 
var endDate16 = (year16)+'-12-31'; 

var l816 = L8.filterDate(startDate16,endDate16)
    .map(maskL8sr)
    .map(addIndicesL8)

var composite16 = l816
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask16 = composite16.select('MNDWI').lt(0.07);
var compositeNew16 = composite16
                        .updateMask(MNDWIMask16)
                        .updateMask(elevationMask);
                        
// Map.addLayer(compositeNew16.clip(East_Luwu), visParL8, 'Landsat 8 Composite 2016')


///APLLIED MODEL to 2016
var image16 = compositeNew16.select(bands).clip(East_Luwu);
var classifiedrf16 = image16.select(bands).classify(classifier_8); 
                      
///pixel count to reduc noise
var pixelcount16 = classifiedrf16.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask16 = pixelcount16.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask16 = classifiedrf16.select('classification').gt(0);
var classed16= classifiedrf16.updateMask(countmask16).updateMask(classMask16);
// Map.addLayer (classed16, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2016');
var imageex16 = composite16.select('B5','B6','B4');

////////////////////////////////////////
///////////Applied 2017//////////////////
////////////////////////////////////////

var year17 = 2017;
var startDate17 = (year17)+'-01-01'; 
var endDate17 = (year17)+'-12-31'; 

var l817 = L8.filterDate(startDate17,endDate17)
    .map(maskL8sr)
    .map(addIndicesL8);

var composite17 = l817
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask17 = composite17.select('MNDWI').lt(0.07);
var compositeNew17 = composite17
                        .updateMask(MNDWIMask17)
                        .updateMask(elevationMask);
                        
// Map.addLayer(compositeNew16.clip(East_Luwu), visParL8, 'Landsat 8 Composite 2017')


///APLLIED MODEL to 2017
var image17 = compositeNew17.select(bands).clip(East_Luwu);
var classifiedrf17 = image17.select(bands).classify(classifier_8); 
                      
///pixel count to reduc noise
var pixelcount17 = classifiedrf17.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask17 = pixelcount17.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask17 = classifiedrf17.select('classification').gt(0);
var classed17= classifiedrf17.updateMask(countmask17).updateMask(classMask17);
// Map.addLayer (classed17, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2017');
var imageex17 = composite17.select('B5','B6','B4');

////////////////////////////////////////
///////////Applied 2018//////////////////
////////////////////////////////////////

var year18 = 2018;
var startDate18 = (year18)+'-01-01'; 
var endDate18 = (year18)+'-12-31'; 

var l818 = L8.filterDate(startDate18,endDate18)
    .map(maskL8sr)
    .map(addIndicesL8);

var composite18 = l818
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask18 = composite18.select('MNDWI').lt(0.07);
var compositeNew18 = composite18
                        .updateMask(MNDWIMask18)
                        .updateMask(elevationMask);
                        
// Map.addLayer(compositeNew18.clip(East_Luwu), visParL8, 'Landsat 8 Composite 2018')


///APLLIED MODEL to 2018
var image18 = compositeNew18.select(bands).clip(East_Luwu);
var classifiedrf18 = image18.select(bands).classify(classifier_8); 
                      
///pixel count to reduc noise
var pixelcount18 = classifiedrf18.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask18 = pixelcount18.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask18 = classifiedrf18.select('classification').gt(0);
var classed18= classifiedrf18.updateMask(countmask18).updateMask(classMask18);
// Map.addLayer (classed18, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2018');
var imageex18 = composite18.select('B5','B6','B4');

////////////////////////////////////////
///////////Applied 2019//////////////////
////////////////////////////////////////

var year19 = 2019;
var startDate19 = (year19)+'-01-01'; 
var endDate19 = (year19)+'-12-31'; 

var l819 = L8.filterDate(startDate19,endDate19)
    .map(maskL8sr)
    .map(addIndicesL8);

var composite19 = l819
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask19 = composite19.select('MNDWI').lt(0.07);
var compositeNew19 = composite19
                        .updateMask(MNDWIMask19)
                        .updateMask(elevationMask);
                        
// Map.addLayer(compositeNew19.clip(East_Luwu), visParL8, 'Landsat 8 Composite 2019')


///APLLIED MODEL to 2019
var image19 = compositeNew19.select(bands).clip(East_Luwu);
var classifiedrf19 = image19.select(bands).classify(classifier_8); 
                      
///pixel count to reduc noise
var pixelcount19 = classifiedrf19.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask19 = pixelcount19.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask19 = classifiedrf19.select('classification').gt(0);
var classed19= classifiedrf19.updateMask(countmask19).updateMask(classMask19);
// Map.addLayer (classed19, {min: 1, max: 1, palette:'green'}, 'Mangrove Extent 2019');
var imageex19 = composite19.select('B5','B6','B4');

////////////////////////////////////////
///////////Applied 2020//////////////////
////////////////////////////////////////

var year20 = 2020;
var startDate20 = (year20)+'-01-01'; 
var endDate20 = (year20)+'-12-31'; 

var l820 = L8.filterDate(startDate20,endDate20)
    .map(maskL8sr)
    .map(addIndicesL8);

var composite20 = l820
              .median() //or ,qualityMosaic('NDVI')
              .clip(East_Luwu); 

///CLIP SRTM and MNDWI
var MNDWIMask20 = composite20.select('MNDWI').lt(0.07);
var compositeNew20 = composite20
                        .updateMask(MNDWIMask20)
                        .updateMask(elevationMask);
                        
Map.addLayer(compositeNew20.clip(East_Luwu), visParL8, 'Landsat 8 Composite 2020');


///APLLIED MODEL to 2020
var image20 = compositeNew20.select(bands).clip(East_Luwu);
var classifiedrf20 = image20.select(bands).classify(classifier_8); 
                      
///pixel count to reduc noise
var pixelcount20 = classifiedrf20.connectedPixelCount(100, false); //Create an image that shows the number of pixels each pixel is connected to
var countmask20 = pixelcount20.select(0).gt(6.25); //filter out all pixels connected to 4 or less 

//Mask the results to only display mangrove extent
var classMask20 = classifiedrf20.select('classification').gt(0);
var classed20= classifiedrf20.updateMask(countmask20).updateMask(classMask20).clip(coastal);
Map.addLayer (classed20, {min: 1, max: 1, palette:'pink'}, 'Mangrove Extent 2020');
var imageex20 = composite20.select('B5','B6','B4');

/////////////////////////////// -- End of Landsat-8 Processing -- ///////////////////////////////


////////////-----------------------------------------------------//////////
//////////// --   Mangrove area calculation for 2000 - 2020   -- //////////
////////////-----------------------------------------------------//////////

//Use reduceRegion with a Sum reducer to calculate total area
var get2000 = classed00.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2001 = classed01.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');
      
var get2002 = classed02.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2003 = classed03.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');
      
var get2004 = classed04.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');
      
var get2005 = classed05.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2006 = classed06.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2007 = classed07.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');
      
var get2008 = classed08.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');
      
var get2009 = classed09.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2010 = classed10.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2011 = classed11.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');
      
var get2012 = classed12.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2013 = classed13.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2014 = classed14.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');
      
var get2015= classed15.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');
      
var get2016 = classed16.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2017 = classed17.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2018 = classed18.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2019 = classed19.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');

var get2020 = classed20.multiply(ee.Image.pixelArea()).divide(10000).reduceRegion({
      reducer:ee.Reducer.sum(),
      geometry:coastal,
      scale: 30,
      maxPixels:1e13,
      tileScale: 16
      }).get('classification');
      
print(get2000, 'Mangrove Extent 2000 in ha');
print(get2001, 'Mangrove Extent 2001 in ha');
print(get2002, 'Mangrove Extent 2002 in ha');
print(get2003, 'Mangrove Extent 2003 in ha');
print(get2004, 'Mangrove Extent 2004 in ha');
print(get2005, 'Mangrove Extent 2005 in ha');
print(get2006, 'Mangrove Extent 2006 in ha');
print(get2007, 'Mangrove Extent 2007 in ha');
print(get2008, 'Mangrove Extent 2008 in ha');
print(get2009, 'Mangrove Extent 2009 in ha');
print(get2010, 'Mangrove Extent 2010 in ha');
print(get2011, 'Mangrove Extent 2011 in ha');
print(get2012, 'Mangrove Extent 2012 in ha');
print(get2013, 'Mangrove Extent 2013 in ha');
print(get2014, 'Mangrove Extent 2014 in ha');
print(get2015, 'Mangrove Extent 2015 in ha');
print(get2016, 'Mangrove Extent 2016 in ha');
print(get2017, 'Mangrove Extent 2017 in ha');
print(get2018, 'Mangrove Extent 2018 in ha');
print(get2019, 'Mangrove Extent 2019 in ha');
print(get2020, 'Mangrove Extent 2020 in ha');


// //Export Mangrove
// Export.image.toDrive({
//   image: classed00,
//   description: 'Mangrove_Extent_2000_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed01,
//   description: 'Mangrove_Extent_2001_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed02,
//   description: 'Mangrove_Extent_2002_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed03,
//   description: 'Mangrove_Extent_2003_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed04,
//   description: 'Mangrove_Extent_2004_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed05,
//   description: 'Mangrove_Extent_2005_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed06,
//   description: 'Mangrove_Extent_2006_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed07,
//   description: 'Mangrove_Extent_2007_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed08,
//   description: 'Mangrove_Extent_2008_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed09,
//   description: 'Mangrove_Extent_2009_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed10,
//   description: 'Mangrove_Extent_2010_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed11,
//   description: 'Mangrove_Extent_2011_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed12,
//   description: 'Mangrove_Extent_2012_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed13,
//   description: 'Mangrove_Extent_2013_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed14,
//   description: 'Mangrove_Extent_2014_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed15,
//   description: 'Mangrove_Extent_2015_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed16,
//   description: 'Mangrove_Extent_2016_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed17,
//   description: 'Mangrove_Extent_2017_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed18,
//   description: 'Mangrove_Extent_2018_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: classed19,
//   description: 'Mangrove_Extent_2019_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });
  
// Export.image.toDrive({
//   image: classed20,
//   description: 'Mangrove_Extent_2020_Filter',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// //Export image
// Export.image.toDrive({
//   image: imageex00,
//   description: 'imageex00',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex01,
//   description: 'imageex01',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex02,
//   description: 'imageex02',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex03,
//   description: 'imageex03',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex04,
//   description: 'imageex04',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex05,
//   description: 'imageex05',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex06,
//   description: 'imageex06',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex07,
//   description: 'imageex07',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex08,
//   description: 'imageex08',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex09,
//   description: 'imageex09',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex10,
//   description: 'imageex10',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex11,
//   description: 'imageex11',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex12,
//   description: 'imageex12',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex13,
//   description: 'imageex13',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex14,
//   description: 'imageex14',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex15,
//   description: 'imageex15',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex16,
//   description: 'imageex16',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex17,
//   description: 'imageex17',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex18,
//   description: 'imageex18',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// Export.image.toDrive({
//   image: imageex19,
//   description: 'imageex19',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });
  
// Export.image.toDrive({
//   image: imageex20,
//   description: 'imageex20',
//   region: East_Luwu,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });

// var dataset = ee.ImageCollection('LANDSAT/MANGROVE_FORESTS');

// var mangrovesVis = {
//   min: 0,
//   max: 1.0,
//   palette: ['green'],
// };

// var mangrove = dataset.first().clip(coastal)
// Export.image.toDrive({
//   image: mangrove,
//   description: 'Refmangrove_2000',
//   region: coastal,
//   folder: 'MangroveTS',
//   scale: 30,
//   maxPixels: 1e13
//   });
// Map.addLayer(mangrove, mangrovesVis, 'Mangroves');
