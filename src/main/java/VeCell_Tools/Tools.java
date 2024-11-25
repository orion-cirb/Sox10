package VeCell_Tools;

import VeCell_Tools.Cellpose.CellposeTaskSettings;
import VeCell_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.ImageCalculator;
import ij.plugin.RGBStackMerge;
import ij.plugin.RoiEnlarger;
import ij.plugin.RoiScaler;
import ij.plugin.filter.Analyzer;
import ij.plugin.frame.RoiManager;
import ij.process.AutoThresholder;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.Measure2Distance;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationDistance;
import mcib3d.geom2.measurementsPopulation.PairObjects3DInt;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.spatial.descriptors.G_Function;
import mcib3d.spatial.descriptors.SpatialDescriptor;
import mcib3d.spatial.sampler.SpatialModel;
import mcib3d.spatial.sampler.SpatialRandomHardCore;
import mcib3d.utils.ThreadUtil;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import net.haesleinhuepf.clijx.bonej.BoneJSkeletonize3D;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.Edge;
import sc.fiji.analyzeSkeleton.Graph;
import sc.fiji.analyzeSkeleton.Point;
import sc.fiji.analyzeSkeleton.SkeletonResult;


/**
 * @author ORION-CIRB
 */
public class Tools {
    
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    private final String helpUrl = "https://github.com/orion-cirb/VeCell/tree/version5";
    public CLIJ2 clij2 = CLIJ2.getInstance();
    
    private BufferedWriter resultsGlobal;
    private BufferedWriter resultsDetail;
    
    private Calibration cal;
    private double pixelVol;
    
    String[] chDialog = new String[]{"Astrocytes", "Vessels (optional)"};
    public int imgSeries = 0;
    public int roiScaling = 9;
    public double roiDilation = 50; // um
    
    public String cellposeEnvDir = IJ.isWindows()? System.getProperty("user.home")+File.separator+"miniconda3"+File.separator+"envs"+File.separator+"CellPose" : "/opt/miniconda3/envs/cellpose";
    public final String cellposeModelDir = IJ.isWindows()? System.getProperty("user.home")+File.separator+".cellpose"+File.separator+"models"+File.separator : "";
    public String cellposeModel = "cyto_sox9_naomie";
    public int cellposeDiam = 20;
    public double cellposeStitchTh = 0.5;
    public double minCellVol = 300; // um3
    public double maxCellVol = 3000; // um3
    
    private int nbNei = 10; // K-nearest neighbors
    private boolean computeGFunction = false;
    private int nbRandomSamples = 50;
    
    private double dog1Sigma1 = 4;
    private double dog1Sigma2 = 8;
    private String vesselThMet1 = "Triangle";
    private boolean dog2 = true;
    private double dog2Sigma1 = 7;
    private double dog2Sigma2 = 14;
    private String vesselThMet2 = "Triangle";
    public double minVesselVol = 600; // um3
    private double minVesselLength = 10; // um
    

    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Flush and close an image
     */
    public void closeImage(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom2.Object3DInt");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Get extension of the first image found in the folder
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd" :
                   ext = fileExt;
                   break;
                case "nd2" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }
    
    
    /**
     * Get images with given extension in folder
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + f);
        }
        Collections.sort(images);
        return(images);
    }
       
    
    /**
     * Get image calibration
     */
    public void findImageCalib(IMetadata meta) {
        cal = new Calibration();
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
    }
    
    
    /**
     * Get channels name and add None at the end of channels list
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels(String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs+1];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        channels[chs] = "None";
        return(channels);     
    }
    
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] chMeta) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 180, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channels", new Font("Monospace", Font.BOLD, 12), Color.blue);
        for (int n = 0; n < chDialog.length; n++)
            gd.addChoice(chDialog[n]+": ", chMeta, chMeta[n]);
        
        gd.addMessage("ROIs", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addNumericField("Scaling factor: ", roiScaling, 0);
        
        gd.addMessage("Astrocytes detection", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addStringField("Cellpose model: ", cellposeModel);
        gd.addToSameRow();
        gd.addNumericField("Cellpose diameter: ", cellposeDiam, 0);
        gd.addNumericField("Min volume (µm3): ", minCellVol, 2);
        gd.addToSameRow();
        gd.addNumericField("Max volume (µm3): ", maxCellVol, 2);
        
        gd.addMessage("Astrocytes spatial distribution", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addNumericField("Neighbors nb: ", nbNei, 0);
        gd.addCheckbox("Compute G-function", computeGFunction);
        gd.addToSameRow();
        gd.addNumericField("Random samples nb: ", nbRandomSamples, 0);
        
        String[] methods = AutoThresholder.getMethods();
        gd.addMessage("Vessels segmentation", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addMessage("DoG filter 1", new Font("Monospace", Font.PLAIN, 12), Color.blue);
        gd.addNumericField("Sigma 1: ", dog1Sigma1, 2);
        gd.addToSameRow();
        gd.addNumericField("Sigma 2: ", dog1Sigma2, 2);
        gd.addChoice("Threshold: ", methods, vesselThMet1);
        gd.addMessage("DoG filter 2", new Font("Monospace", Font.PLAIN, 12), Color.blue);
        gd.addCheckbox("Apply filter", dog2);        
        gd.addNumericField("Sigma 1: ", dog2Sigma1, 2);
        gd.addToSameRow();
        gd.addNumericField("Sigma 2: ", dog2Sigma2, 2);
        gd.addChoice("Threshold: ", methods, vesselThMet2);
        gd.addNumericField("Min volume (µm3): ", minVesselVol, 2);
        gd.addToSameRow();
        gd.addNumericField("Min length (µm): ", minVesselLength, 2);
        
        gd.addMessage("Image calibration", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addNumericField("XY pixel size (µm): ", cal.pixelWidth, 4);
        gd.addToSameRow();
        gd.addNumericField("Z pixel size (µm): ", cal.pixelDepth, 4);
        
        gd.addHelp(helpUrl);
        gd.showDialog();
        
        String[] chOrder = new String[chDialog.length];
        for (int n = 0; n < chOrder.length; n++)
            chOrder[n] = gd.getNextChoice();
        
        roiScaling = (int) gd.getNextNumber();
        if (roiScaling <= 0) {
            roiScaling = 1;
            print("WARNING: ROIs cannot be scaled by zero or negative values, ROI scaling factor set to 1");
        }
        
        cellposeModel = gd.getNextString();
        cellposeDiam = (int) gd.getNextNumber();
        minCellVol = gd.getNextNumber();
        maxCellVol = gd.getNextNumber();
        
        nbNei = (int)gd.getNextNumber();
        computeGFunction = gd.getNextBoolean();
        nbRandomSamples = (int)gd.getNextNumber();
        
        dog1Sigma1 = (int) gd.getNextNumber();
        dog1Sigma2 = (int) gd.getNextNumber();
        vesselThMet1 = gd.getNextChoice();
        dog2 = gd.getNextBoolean();
        dog2Sigma1 = (int) gd.getNextNumber();
        dog2Sigma2 = (int) gd.getNextNumber();
        vesselThMet2 = gd.getNextChoice();
        minVesselVol = gd.getNextNumber();
        minVesselLength = gd.getNextNumber();
        
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        pixelVol = cal.pixelWidth * cal.pixelHeight * cal.pixelDepth;
        
        if (gd.wasCanceled())
            chOrder = null;
        return(chOrder);
    }
    
    
    /**
     * Save images specific channel before sending it to QuantileBasedNormalization plugin
     * @throws Exception
     */
    public void saveChannel(ArrayList<String> imgFiles, int series, int channel, String normDir, String extension) throws Exception {
        try {           
            for (String f : imgFiles) {
                String imgName = FilenameUtils.getBaseName(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setQuiet(true);
                
                // Open and save vessels channel
                ImagePlus imgVessels = openChannel(options, series, channel);
                IJ.saveAs(imgVessels, "Tiff", normDir+imgName+extension);
                closeImage(imgVessels);
            }
        } catch (Exception e) {
            throw e; 
        }
    }
    
    
    public ImagePlus openChannel(ImporterOptions options, int series, int channel) throws FormatException, IOException {
        options.setCBegin(series, channel);
        options.setCEnd(series, channel);
        ImagePlus img = BF.openImagePlus(options)[0];
        return(img);
    }
    
    
    /**
     * Delete images specific channel after it was sent to QuantileBasedNormalization plugin
     * @throws Exception
     */
    public void deleteChannel(String dir, ArrayList<String> imageFiles, String extension) throws Exception {
        try {
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                new File(dir+rootName+extension).delete();
            }
        } catch (Exception e) {
            throw e; 
        }
    }
    
    
    public List<Roi> loadRois(String imageDir, String rootName) {
        String roiName = imageDir+rootName;
        roiName = new File(roiName+".zip").exists() ? roiName+".zip" : roiName+".roi";
        if (new File(roiName).exists()) {
            RoiManager rm = new RoiManager(false);
            rm.runCommand("Open", roiName);
            List<Roi> rois = Arrays.asList(rm.getRoisAsArray());
            return(rois);
        } else {
            print("ERROR: No ROI file found for image " + rootName + ", image not analyzed");
            return(null);
        }
    }
        

    /**
     * Scale ROIs
     */
    public List<Roi> scaleRois(List<Roi> rois, int scale) {
        List<Roi> scaledRois = new ArrayList<Roi>();
        for (Roi roi : rois) {
            Roi scaledRoi = new RoiScaler().scale(roi, scale, scale, false);
            Rectangle rect = roi.getBounds();
            scaledRoi.setLocation(rect.x*scale, rect.y*scale);
            scaledRoi.setName(roi.getName());
            scaledRois.add(scaledRoi);
        }
        return scaledRois;
    }
    
    
    /**
     * Get bounding box of multiple dilated ROIs combined together
     */
    public Roi getBoundingBox(List<Roi> rois) {
        Rectangle rect0 = RoiEnlarger.enlarge(rois.get(0), roiDilation/cal.pixelWidth).getBounds();
        int minX = rect0.x, minY = rect0.y, maxX = rect0.x+rect0.width, maxY = rect0.y+rect0.height;
        for (Roi roi : rois) {
            Rectangle rect = RoiEnlarger.enlarge(roi, roiDilation/cal.pixelWidth).getBounds();
            if(rect.x < minX) minX = rect.x;
            if(rect.y < minY) minY = rect.y;
            if(rect.x+rect.width > maxX) maxX = rect.x+rect.width;
            if(rect.y+rect.height > maxY) maxY = rect.y+rect.height;
        }
        minX = Math.max(minX, 0);
        minY = Math.max(minY, 0);
        return new Roi(minX, minY, maxX-minX, maxY-minY);
    }
    
    
    /**
     * Translate ROIs
     */
    public void translateRois(List<Roi> rois, Roi bBox) {
        Rectangle rectBBox = bBox.getBounds();
        for (Roi roi: rois) {
            Rectangle rectRoi = roi.getBounds();
            roi.setLocation(rectRoi.x-rectBBox.x, rectRoi.y-rectBBox.y);
        }
    }
    
    
    /**
     * Save ROIs in a zip file
     */
    public void saveRois(List<Roi> rois, String outDir, String imgName) {
        RoiManager rm = new RoiManager(false);
        for(Roi roi: rois) {
            String layerName = roi.getName().split(" ")[1];
            rm.addRoi(roi);
            String roiName = rm.getName(rm.getCount()-1) + " " + layerName;
            rm.rename(rm.getCount()-1, roiName);
        }
        rm.deselect();
        rm.save(outDir+imgName+"_rois.zip");
        rm.close();
    }
    
    
    /**
     * Crop stack around ROI
     */
    public ImagePlus cropImage(ImagePlus img, Roi roi) {
        img.setRoi(roi);
        ImagePlus imgCrop = img.crop("stack");
        img.deleteRoi();
        return(imgCrop);
    }
    
    
    /*
     * Look for all 3D cells in a Z-stack: 
     * - apply Cellpose 2D slice by slice 
     * - let CellPose reconstruct cells in 3D using the stitch_threshold parameter
     */
    public ImagePlus cellposeDetection(ImagePlus imgIn) throws IOException{
        ImagePlus img = imgIn.duplicate();
       
        // Define Cellpose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(cellposeModel.equals("cyto")? cellposeModel : cellposeModelDir+cellposeModel, 1, cellposeDiam, cellposeEnvDir); // need to add Cellpose models folder path if own model (for Windows only, not Linux)
        settings.setStitchThreshold(cellposeStitchTh);
        settings.useGpu(true);
       
        // Run Cellpose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, img);
        ImagePlus imgOut = cellpose.run();
        if(imgOut.getNChannels() > 1) 
            imgOut.setDimensions(1, imgOut.getNChannels(), 1);
        imgOut.setCalibration(cal);
       
        closeImage(img);
        return imgOut;
    }
    
    
    /**
     * Detect vessels applying a median filter + DoG filter + threshold + closing filter + median filter
     */
    public ImagePlus vesselsSegmentation(ImagePlus img) {
        ImagePlus imgMed = medianFilter3D(img, 4, 1);
        ImagePlus imgDOG = dogFilter2D(imgMed, dog1Sigma1, dog1Sigma2);
        ImagePlus imgBin = threshold(imgDOG, vesselThMet1);
        if(dog2) {
            ImagePlus imgDOG2 = dogFilter2D(imgMed, dog2Sigma1, dog2Sigma2);
            ImagePlus imgBin2 = threshold(imgDOG2, vesselThMet2);
            new ImageCalculator().run("Max stack", imgBin, imgBin2);
            closeImage(imgDOG2);
            closeImage(imgBin2);
        }
        ImagePlus imgClose = closingFilter3D(imgBin, 8, 1);
        ImagePlus imgOut = medianFilter3D(imgClose, 1, 1);

        ImageInt imgLabels = new ImageLabeller().getLabels(ImageHandler.wrap(imgOut));
        imgLabels.setCalibration(cal);
        
        Objects3DIntPopulation pop = new Objects3DIntPopulation(imgLabels);
        popFilterOneZ(pop);
        System.out.println(pop.getNbObjects() + " vessels detected");
        popFilterVol(pop, minVesselVol, Double.MAX_VALUE);
        System.out.println(pop.getNbObjects() + " vessels remaining after size filtering");
        
        ImageHandler imgFilter = ImageHandler.wrap(img).createSameDimensions();
        for(Object3DInt obj: pop.getObjects3DInt())
            obj.drawObject(imgFilter, 255);
        imgFilter.setCalibration(cal);
        
        closeImage(imgMed);
        closeImage(imgDOG);
        closeImage(imgBin);
        closeImage(imgClose);
        closeImage(imgOut);
        closeImage(imgLabels.getImagePlus());
        
        return(imgFilter.getImagePlus());
    }
          
    
    /**
     * Median filter using CLIJ
     */ 
    public ImagePlus medianFilter3D(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median3DSphere(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
       clij2.release(imgCL);
       ImagePlus imgMed = clij2.pull(imgCLMed);
       clij2.release(imgCLMed);
       return(imgMed);
    }
    
      
    /**
     * Difference of Gaussians using CLIJ
     */ 
    public ImagePlus dogFilter2D(ImagePlus img, double size1, double size2) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian2D(imgCL, imgCLDOG, size1, size1, size2, size2);
        clij2.release(imgCL);
        ImagePlus imgDOG = clij2.pull(imgCLDOG); 
        clij2.release(imgCLDOG);
        return(imgDOG);
    }
    
    
    /**
     * Automatic thresholding using CLIJ2
     */
    private ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCL);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
    
    /**
     * Closing filtering using CLIJ2
     */ 
    private ImagePlus closingFilter3D(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMax = clij2.create(imgCL);
       clij2.maximum3DSphere(imgCL, imgCLMax, sizeXY, sizeXY, sizeZ);

       ClearCLBuffer imgCLMin = clij2.create(imgCLMax);
       clij2.minimum3DSphere(imgCLMax, imgCLMin, sizeXY, sizeXY, sizeZ);
       ImagePlus imgMin = clij2.pull(imgCLMin);
       
       clij2.release(imgCL);
       clij2.release(imgCLMax);
       clij2.release(imgCLMin);
       return(imgMin);
    }
    
    
    /**
     * Compute (inverse) 3D distance map
     */
    public ImageFloat distanceMap3D(ImagePlus img, boolean inverse) {
        img.setCalibration(cal);
        ImageFloat edt = new EDT().run(ImageHandler.wrap(img), 0, inverse, ThreadUtil.getNbCpus());
        return(edt);
    }
    
    
    /**
     * Skeletonize 3D with CLIJ2
     */
    public ImagePlus skeletonize3D(ImagePlus img) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLSkel = clij2.create(imgCL);
        new BoneJSkeletonize3D().bonejSkeletonize3D(clij2, imgCL, imgCLSkel);
        ImagePlus imgSkel = clij2.pull(imgCLSkel);
        clij2.release(imgCL);
        clij2.release(imgCLSkel);
        
        IJ.run(imgSkel, "8-bit","");
        imgSkel.setCalibration(cal);
        return(imgSkel);
    }
    
    
    /**
     * Prune skeleton branches with length smaller than threshold
     * https://imagej.net/plugins/analyze-skeleton/
     */
    public ImagePlus pruneSkeleton(ImagePlus image) {
        // Analyze skeleton
        AnalyzeSkeleton_ skel = new AnalyzeSkeleton_();
        skel.setup("", image);
        SkeletonResult skelResult = skel.run(AnalyzeSkeleton_.NONE, false, false, null, true, false);

        // Create copy of input image
        ImagePlus prunedImage = image.duplicate();
        if(skelResult.getBranches() != null) {
            ImageStack outStack = prunedImage.getStack();

            // Get graphs (one per skeleton in the image)
            Graph[] graphs = skelResult.getGraph();

            // Get list of end-points
            ArrayList<Point> endPoints = skelResult.getListOfEndPoints();

            for(Graph graph: graphs) {
                ArrayList<Edge> listEdges = graph.getEdges();

                // Go through all branches and remove branches under threshold in duplicate image
                for(Edge e: listEdges) {
                    ArrayList<Point> p1 = e.getV1().getPoints();
                    boolean v1End = endPoints.contains(p1.get(0));
                    ArrayList<Point> p2 = e.getV2().getPoints();
                    boolean v2End = endPoints.contains(p2.get(0));
                    // If any of the vertices is end-point 
                    if(v1End || v2End) {
                        if(e.getLength() < minVesselLength) { // in microns
                            if(v1End)
                                outStack.setVoxel(p1.get(0).x, p1.get(0).y, p1.get(0).z, 0);
                            if(v2End)
                                outStack.setVoxel(p2.get(0).x, p2.get(0).y, p2.get(0).z, 0);
                            for(Point p: e.getSlabs())
                                outStack.setVoxel(p.x, p.y, p.z, 0);
                        }
                    }
                }
            }
        }
        
        return(prunedImage);
    }
    
    
    /**
     * Get object in (dilated) ROI
     * @param dilationFactor in microns
     */
    public ImagePlus getImgInRoi(ImagePlus img, Roi roi, double dilationFactor) {
        ImagePlus imgClear = img.duplicate();
        imgClear.setRoi(dilationFactor == 0? roi : RoiEnlarger.enlarge(roi, dilationFactor/cal.pixelWidth)); // pixels
        IJ.run(imgClear, "Clear Outside", "stack");
        imgClear.setCalibration(cal);
        return(imgClear);
    }
    
    
    /**
     * Get object in (dilated) ROI
     */
    public Object3DInt getObjInRoi(ImagePlus img, Roi roi, double dilationFactor) {
        ImagePlus imgClear = getImgInRoi(img, roi, dilationFactor);
        Object3DInt obj = new Object3DInt(ImageHandler.wrap(imgClear));
        closeImage(imgClear);
        return(obj);
    }
    
    
    /**
     * Return population in ROI
     */
    public Objects3DIntPopulation getPopInRoi(ImagePlus img, Roi roi) throws IOException{
        ImagePlus imgClear = getImgInRoi(img, roi, 20);

        // Filter detections
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageInt.wrap(imgClear));
        popFilterOneZ(pop);
        popFilterCentroid(pop, roi);
        System.out.println(pop.getNbObjects() + " cells detected in ROI");
        popFilterVol(pop, minCellVol, maxCellVol);
        System.out.println(pop.getNbObjects() + " cells remaining after size filtering");
        pop.resetLabels();

        closeImage(imgClear);
        return(pop);
    }
     
    
    /**
     * Remove objects that appear in only one z-slice
     */
    private void popFilterOneZ(Objects3DIntPopulation pop) {
        pop.getObjects3DInt().removeIf(p -> (p.getObject3DPlanes().size() == 1));
        pop.resetLabels();
    }
    
    
    /**
     * Remove objects that do not have their centroid into a given ROI
     */
    public void popFilterCentroid(Objects3DIntPopulation pop, Roi roi) {
        pop.getObjects3DInt().removeIf(p -> (!roi.contains(new MeasureCentroid(p).getCentroidAsPoint().getRoundX(), 
                                                           new MeasureCentroid(p).getCentroidAsPoint().getRoundY())));
        pop.resetLabels();
    }
    
    
    /**
     * Remove objects with volume < min and volume > max
     */
    private void popFilterVol(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
    
    /**
     * Remove objects with intensity < min
     */
    public void popFilterInt(Objects3DIntPopulation pop, ImagePlus img, double minInt) {
        ImageHandler imh = ImageHandler.wrap(img);
        pop.getObjects3DInt().removeIf(p -> (new MeasureIntensity(p, imh).getValueMeasurement(MeasureIntensity.INTENSITY_AVG) < minInt));
        pop.resetLabels();
    }
    
    
    /**
     * Compute distance between each cell and its closest vessel
     */
    public void computeCellsVesselsDists(Objects3DIntPopulation cellPop, ImageFloat vesselDistMapInv) {
        for (Object3DInt cell: cellPop.getObjects3DInt()) {
            double dist = vesselDistMapInv.getPixel(new MeasureCentroid​(cell).getCentroidAsPoint());
            cell.setCompareValue(dist);
        }
    }
    
    
    /**
     * Draw results obtained in ROI
     */
    public void drawResultsInRoi(Objects3DIntPopulation popCell, Object3DInt objVessel, Object3DInt objSkel, ImageHandler imhCell,
                                 ImageHandler imhCellDist, ImageHandler imhVessel, ImageHandler imhSkel) {
        // Cells
        popCell.drawInImage(imhCell);
        
        // Cells + vessels
        if (objVessel != null) {
            for (Object3DInt cell: popCell.getObjects3DInt()) {
                cell.drawObject(imhCellDist, (float)cell.getCompareValue()+10);
            }
            objVessel.drawObject(imhVessel, 255);
            objSkel.drawObject(imhSkel, 255);
        }
    }
    
    
    /**
     * Save results in images
     */
    public void saveCloseDrawings(ImagePlus imgCell, ImagePlus imgVessel, ImageHandler imhCell, ImageHandler imhCellDist, 
                            ImageHandler imhVessel, ImageHandler imhSkel, String outDir, String imgName) {
        // Cells
        ImagePlus[] imgColors1 = {imhCell.getImagePlus(), null, null, imgCell};
        ImagePlus imgObjects1 = new RGBStackMerge().mergeHyperstacks(imgColors1, true);
        imgObjects1.setCalibration(cal);
        new FileSaver(imgObjects1).saveAsTiff(outDir+imgName+"_cells.tif"); 
        closeImage(imgObjects1);
        
        // Cells + vessels
        if (imgVessel != null) {
            IJ.run(imhCellDist.getImagePlus(), "mpl-inferno", "");
            IJ.run(imhSkel.getImagePlus(), "Green", "");
            IJ.run(imhVessel.getImagePlus(), "Blue", "");
            
            ImagePlus[] imgColors2 = {imhCellDist.getImagePlus(), imhSkel.getImagePlus(), imhVessel.getImagePlus(), imgVessel};
            ImagePlus imgObjects2 = new RGBStackMerge().mergeHyperstacks(imgColors2, true);
            imgObjects2.setCalibration(cal);
            new FileSaver(imgObjects2).saveAsTiff(outDir+imgName+"_vessels.tif");
            closeImage(imgObjects2);
        }
        
        closeImage(imhCell.getImagePlus());
        closeImage(imhCellDist.getImagePlus());
        closeImage(imhVessel.getImagePlus());
        closeImage(imhSkel.getImagePlus());
    }
    
    
    /**
     * Write headers in results files
     */
    public void writeHeaders(String outDir, boolean computeVessel) throws IOException {
        // Global results
        FileWriter fileGlobal = new FileWriter(outDir + "globalResults.csv", false);
        resultsGlobal = new BufferedWriter(fileGlobal);
        resultsGlobal.write("Image name\tROI name\tROI area (µm²)\tROI volume (µm³)\tNb cells\tCells total volume (µm³)\t"
                + "Cells mean volume (µm³)\tCells volume SD (µm³)\tCells mean distance to closest neighbor (µm)\t"
                + "Cells distance to closest neighbor SD (µm)\tCells mean distance to "+nbNei+" closest neighbors (µm)\t"
                + "Cells distance to "+nbNei+" closest neighbors SD (µm)"+"\tCells mean of max distance to "+nbNei+" neighbors (µm)\t"
                + "Cells SD of max distance to "+nbNei+" neighbors (µm)");
        if (computeGFunction) resultsGlobal.write("\tCells G-function SDI\tCells G-function AUC difference");
        if (computeVessel) resultsGlobal.write("\tCells mean distance to closest vessel (µm)\tVessels total volume (µm³)"
                + "\tVessels total length (µm)\tBranches mean length (µm)\tNb branches\tNb junctions\tVessels mean diameter (µm)\tVessels diameter SD (µm)");
        resultsGlobal.write("\n");
        resultsGlobal.flush();

        // Detailed results
        FileWriter fileDetail = new FileWriter(outDir +"cellsResults.csv", false);
        resultsDetail = new BufferedWriter(fileDetail);
        resultsDetail.write("Image name\tROI name\tCell ID\tCell volume (µm³)\tCell distance to closest neighbor (µm)\t"
                + "Cell mean distance to "+nbNei+" closest neighbors (µm)\tCell max distance to "+nbNei+" closest neighbors (µm)");
        if(computeVessel) resultsDetail.write("\tCell distance to closest vessel (µm)\tClosest vessel diameter (µm)");
        resultsDetail.write("\n");
        resultsDetail.flush();
    }
    
    
    /**
     * Compute parameters and save them in results files
     * @throws java.io.IOException
     */
    public void writeResultsInRoi(Objects3DIntPopulation popCellRoi, Object3DInt objSkelRoiDil, Object3DInt objSkelRoi,
                                  ImageFloat distMap, ImagePlus imgCell, Roi roi, double vesselVol, String outDir, String imgName) throws IOException {
        
        DescriptiveStatistics cellsVolume = new DescriptiveStatistics();
        DescriptiveStatistics cellsClosestVesselDist = new DescriptiveStatistics();
        DescriptiveStatistics cellsClosestNeighborDist = new DescriptiveStatistics();
        DescriptiveStatistics cellsNeighborsMeanDist = new DescriptiveStatistics();
        DescriptiveStatistics cellsNeighborsMaxDist = new DescriptiveStatistics();
        
        // CELLS INDIVIDUAL STATISTICS
        print("Computing cells individual statistics...");
        MeasurePopulationDistance allCellsDists = new MeasurePopulationDistance​(popCellRoi, popCellRoi, Double.POSITIVE_INFINITY, "DistCenterCenterUnit");
        for (Object3DInt cell: popCellRoi.getObjects3DInt()) {
            double cellVol = new MeasureVolume(cell).getVolumeUnit();
            cellsVolume.addValue(cellVol);            
            resultsDetail.write(imgName+"\t"+roi.getName()+"\t"+cell.getLabel()+"\t"+cellVol);
            resultsDetail.flush();
            
            if(popCellRoi.getNbObjects() == 1) {
                resultsDetail.write("\t"+Double.NaN+"\t"+Double.NaN+"\t"+Double.NaN);
                resultsDetail.flush();
            } else {
                List<PairObjects3DInt> cellCellsDists = allCellsDists.getPairsObject1(cell.getLabel(), true);
                double closestNeighborDist = cellCellsDists.get(1).getPairValue();
                cellsClosestNeighborDist.addValue(closestNeighborDist);

                DescriptiveStatistics cellNeighborsDists = new DescriptiveStatistics();
                for (int d=1; d <= Math.min(nbNei, popCellRoi.getNbObjects()-1); d++)
                    cellNeighborsDists.addValue(cellCellsDists.get(d).getPairValue());

                double closestNeighborsMeanDist = cellNeighborsDists.getMean();
                cellsNeighborsMeanDist.addValue(closestNeighborsMeanDist);
                double closestNeighborsMaxDist = cellNeighborsDists.getMax();
                cellsNeighborsMaxDist.addValue(closestNeighborsMaxDist);

                resultsDetail.write("\t"+closestNeighborDist+"\t"+closestNeighborsMeanDist+"\t"+closestNeighborsMaxDist);
                resultsDetail.flush();
            }
            
            cellsClosestVesselDist.addValue(cell.getCompareValue());
            if (objSkelRoiDil != null) {
                double diam = 2*distMap.getPixel(new Measure2Distance(cell, objSkelRoiDil).getBorder2Pix());
                resultsDetail.write("\t"+cell.getCompareValue()+"\t"+diam);
                resultsDetail.flush();
            }
            
            resultsDetail.write("\n");
            resultsDetail.flush();
        }
          
        // CELLS GLOBAL STATISTICS
        print("Computing cells global statistics...");
        double[] roiParams = roiParams(roi, imgCell);        
        resultsGlobal.write(imgName+"\t"+roi.getName()+"\t"+roiParams[0]+"\t"+roiParams[1]+"\t"+popCellRoi.getNbObjects()+"\t"+
                            cellsVolume.getSum()+"\t"+cellsVolume.getMean()+"\t"+cellsVolume.getStandardDeviation()+"\t"+
                            cellsClosestNeighborDist.getMean()+"\t"+cellsClosestNeighborDist.getStandardDeviation()+"\t"+
                            cellsNeighborsMeanDist.getMean()+"\t"+cellsNeighborsMeanDist.getStandardDeviation()+"\t"+
                            cellsNeighborsMaxDist.getMean()+"\t"+cellsNeighborsMaxDist.getStandardDeviation());
        
        if (computeGFunction) {
            if(popCellRoi.getNbObjects() <= 1) {
                resultsGlobal.write("\t"+Double.NaN+"\t"+Double.NaN);
                resultsGlobal.flush();
            } else {
                System.out.println("Computing G-function-related spatial distribution index...");
                Object3DInt mask = roiMask(imgCell, roi);
                double minDist = Math.pow(3*minCellVol/(4*Math.PI*pixelVol), 1/3) * 2; // min distance = 2 * min cell radius (in pixels)
                String plotName = outDir + imgName + "_" + roi.getName() + "_gfunction.tif";
                double[] res = computeSdiG(popCellRoi, mask, imgCell, minDist, nbRandomSamples, plotName);
                resultsGlobal.write("\t"+res[0]+"\t"+res[1]);
                resultsGlobal.flush();
            }
        }
        
        // VESSELS STATISTICS
        if(objSkelRoi == null) {
            resultsGlobal.write("\t0\t0\t0\t0\t0\t0\t0\t0");
            resultsGlobal.flush();
        } else {
            print("Computing vessels statistics...");
            ImageHandler imhSkel = ImageHandler.wrap(imgCell).createSameDimensions();
            objSkelRoi.drawObject(imhSkel, 255);
            IJ.run(imhSkel.getImagePlus(), "8-bit","");

            AnalyzeSkeleton_ analyzeSkeleton = new AnalyzeSkeleton_();
            analyzeSkeleton.setup("", imhSkel.getImagePlus());
            SkeletonResult skelResult = analyzeSkeleton.run(AnalyzeSkeleton_.NONE, false, false, null, true, false);
            closeImage(imhSkel.getImagePlus());

            if(skelResult.getBranches() == null) {
                resultsGlobal.write("\t0\t0\t0\t0\t0\t0\t0\t0");
                resultsGlobal.flush();
            } else {   
                double[] branchLengths = skelResult.getAverageBranchLength();
                int[] branchNumbers = skelResult.getBranches();
                double totalLength = 0;
                for (int i = 0; i < branchNumbers.length; i++)
                    totalLength += branchNumbers[i] * branchLengths[i];

                DescriptiveStatistics diams = new DescriptiveStatistics();
                for (Point pt: skelResult.getListOfSlabVoxels())
                    diams.addValue(2*distMap.getPixel(pt.x, pt.y, pt.z));

                resultsGlobal.write("\t"+cellsClosestVesselDist.getMean()+"\t"+vesselVol+"\t"+totalLength+"\t"+StatUtils.mean(branchLengths)+"\t"+
                                    IntStream.of(branchNumbers).sum()+"\t"+IntStream.of(skelResult.getJunctions()).sum()+"\t"+
                                    diams.getMean()+"\t"+diams.getStandardDeviation());
                resultsGlobal.flush();
            }
        }
        
        resultsGlobal.write("\n");
        resultsGlobal.flush();
    }
    
    
    /**
     * Compute ROI area and volume
     */
    public double[] roiParams(Roi roi, ImagePlus img) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        img.setRoi(poly);
        
        ResultsTable rt = new ResultsTable();
        Analyzer analyzer = new Analyzer(img, Analyzer.AREA, rt);
        analyzer.measure();
        double area = rt.getValue("Area", 0);
        double vol = area * img.getNSlices() * cal.pixelDepth;
        double[] params = {area, vol};
        return(params);
    }
    
    
    // Get ROI as a 3D object
    public Object3DInt roiMask(ImagePlus img, Roi roi) {
        ImagePlus imgMask = img.duplicate();
        imgMask.setRoi(roi);
        for (int n = 1; n <= imgMask.getNSlices(); n++) {
            imgMask.setSlice(n);
            IJ.run(imgMask, "Fill", "stack");
            IJ.run(imgMask, "Clear Outside", "stack");
        }
        imgMask.setCalibration(cal);
        
        Object3DInt mask = new Object3DInt​(ImageHandler.wrap(imgMask));
        closeImage(imgMask);
        return(mask);
    }
   

    /**
     * Compute G-function-related Spatial Distribution Index of cells population in a ROI
     * https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000853
     */
    public double[] computeSdiG(Objects3DIntPopulation popInt, Object3DInt roiInt, ImagePlus img, double distHardCore, int numRandomSamples, String plotName) {
        // Convert Object3DInt & Objects3DIntPopulation objects into Object3D & Objects3DPopulation objects
        ImageHandler imhRoi = ImageHandler.wrap(img).createSameDimensions();
        roiInt.drawObject(imhRoi, 1);
        Object3D roi = new Objects3DPopulation(imhRoi).getObject(0);
        closeImage(imhRoi.getImagePlus());
        ImageHandler imhPop = ImageHandler.wrap(img).createSameDimensions();
        popInt.drawInImage(imhPop);
        Objects3DPopulation pop = new Objects3DPopulation(imhPop);
        closeImage(imhPop.getImagePlus());
        
        // Define spatial descriptor and model
        SpatialDescriptor spatialDesc = new G_Function();     
        SpatialModel spatialModel = new SpatialRandomHardCore(pop.getNbObjects(), distHardCore, roi); // average diameter of a cell in pixels
        SpatialStatistics spatialStatistics = new SpatialStatistics(spatialDesc, spatialModel, numRandomSamples, pop); // nb of samples (randomized organizations simulated to compare with the spatial organization of the cells)
        spatialStatistics.setEnvelope(0.05); // 2.5-97.5% envelope error
        spatialStatistics.setVerbose(false);
        double sdiG = spatialStatistics.getSdi();
        double areaG = spatialStatistics.getAUCDifference();
        
        Plot plotG = spatialStatistics.getPlot();
        plotG.draw();
        plotG.addLabel(0.05, 0.1, "SDI = " + String.format("%.5f", sdiG));
        plotG.addLabel(0.05, 0.15, "Area = " + String.format("%.5f", areaG));
        ImagePlus imgPlot = plotG.getImagePlus();
        FileSaver plotSave = new FileSaver(imgPlot);
        plotSave.saveAsTiff(plotName);
        closeImage(imgPlot); 
        double[] results = {sdiG, areaG};
        return(results);
    }

    
    /**
     * Close results files
     */
    public void closeResults() throws IOException {
       resultsGlobal.close();
       resultsDetail.close();
    }
     
}
