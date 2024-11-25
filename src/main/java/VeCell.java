import VeCell_Tools.QuantileBasedNormalization;
import VeCell_Tools.Tools;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.MetadataTools;
import loci.formats.meta.IMetadata;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;


/**
 * @author ORION-CIRB
 */   
public class VeCell implements PlugIn {
    
    private Tools tools = new Tools();
    
    public void run(String arg) {
        try {
            if (!tools.checkInstalledModules()) {
                return;
            }
            
            String imgDir = IJ.getDirectory("Choose images folder");
            if (imgDir == null) {
                return;
            }
            
            // Find extension of first image in input folder
            String fileExt = tools.findImageType(new File(imgDir));
            // Find all images with corresponding extension in folder
            ArrayList<String> imgFiles = tools.findImages(imgDir, fileExt);
            if (imgFiles.isEmpty()) {
                IJ.showMessage("ERROR", "No image found with " + fileExt + " extension in " + imgDir + " folder");
                return;
            }
            
            // Instantiate metadata and reader
            IMetadata meta = MetadataTools.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imgFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);

            // Find channels name
            String[] chMeta = tools.findChannels(imgFiles.get(0), meta, reader);
            
            // Generate dialog box
            String[] chOrder = tools.dialog(chMeta);
            if (chOrder == null) {
                return;
            } else if(chOrder[0].equals("None")) {
                IJ.showMessage("ERROR", "Astrocytes channel not defined");
                return;
            }
            
            // Create output folder
            String outDir = imgDir + File.separator + "Results_" + new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            if (!Files.exists(Paths.get(outDir))) {
                new File(outDir).mkdir();
            }
            
            // Write header in results file
            tools.writeHeaders(outDir, !chOrder[1].equals("None"));
            
            // If asked in dialog box, normalize vessels channel
            String normDir = imgDir + File.separator + "Normalization" + File.separator;
            if (!chOrder[1].equals("None") && !Files.exists(Paths.get(normDir))) {
                tools.print("--- NORMALIZING IMAGES ---");
                // Create output folder for normalized files
                new File(normDir).mkdir();
                // Save images vessels channel
                tools.saveChannel(imgFiles, tools.imgSeries, ArrayUtils.indexOf(chMeta, chOrder[1]), normDir, "-vessels.tif");
                // Normalize images vessels channel
                new QuantileBasedNormalization().run(normDir, imgFiles, "-vessels");
                // Delete images vessels channel
                tools.deleteChannel(normDir, imgFiles, "-vessels.tif");
                tools.print("Normalization done");
            }
            
            IJ.setForegroundColor(255, 255, 255);
            IJ.setBackgroundColor(0, 0, 0);
                 
            for (String f: imgFiles) {
                reader.setId(f);
                String imgName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + imgName + " ---");
                
                tools.print("- Loading ROIs -");
                List<Roi> rawRois = tools.loadRois(imgDir, imgName);
                if (rawRois == null) continue;
                
                List<Roi> rois = tools.scaleRois(rawRois, tools.roiScaling);
                Roi bBox = tools.getBoundingBox(rois);
                tools.translateRois(rois, bBox);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                
                 // Astrocytes channel
                tools.print("- Opening astrocytes channel -");
                ImagePlus imgCell = tools.openChannel(options, tools.imgSeries, ArrayUtils.indexOf(chMeta, chOrder[0]));
                imgCell = tools.cropImage(imgCell, bBox);
                
                tools.print("- Detecting astrocytes -");
                ImagePlus imgCellMask = tools.cellposeDetection(imgCell);
                
                // Vessels channel
                ImagePlus imgVessel = null, imgVesselMask = null, imgVesselSkel = null;
                ImageFloat imgVesselDistMap = null, imgVesselDistMapInv = null;
                if (!chOrder[1].equals("None")){
                    tools.print("- Opening vessels channel -");
                    imgVessel = IJ.openImage(normDir + imgName + "-vessels-normalized.tif");
                    imgVessel = tools.cropImage(imgVessel, bBox);
                    
                    tools.print("- Detecting vessels -");
                    imgVesselMask = tools.vesselsSegmentation(imgVessel);
                    
                    tools.print("- Computing vessels distance maps -");
                    imgVesselDistMap = tools.distanceMap3D(imgVesselMask, false);
                    imgVesselDistMapInv = tools.distanceMap3D(imgVesselMask, true);
                    
                    tools.print("- Computing vessels skeleton -");
                    imgVesselSkel = tools.skeletonize3D(imgVesselMask);
                    // Prune vessels skeleton small branches
                    imgVesselSkel = tools.pruneSkeleton(imgVesselSkel);
                }
                
                // For each roi, open cropped image
                ImageHandler imhCell = ImageHandler.wrap(imgCell).createSameDimensions();
                ImageHandler imhCellDist = imhCell.createSameDimensions();
                ImageHandler imhVessel = imhCell.createSameDimensions();
                ImageHandler imhSkel = imhCell.createSameDimensions();
                for (Roi roi : rois) {
                    tools.print("- Computing parameters and saving results for ROI " + roi.getName() + " -");

                    // Cells channel
                    Objects3DIntPopulation popCellRoi = tools.getPopInRoi(imgCellMask, roi);
                    
                    // Vessels channel
                    double vesselVolRoi = 0;
                    Object3DInt objVesselRoiDil = null, objSkelRoi = null, objSkelRoiDil = null;
                    if (imgVessel != null){
                        vesselVolRoi = new MeasureVolume(tools.getObjInRoi(imgVesselMask, roi, 0)).getVolumeUnit();
                        objVesselRoiDil = tools.getObjInRoi(imgVesselMask, roi, tools.roiDilation);
                        objSkelRoi = tools.getObjInRoi(imgVesselSkel, roi, 0);
                        objSkelRoiDil = tools.getObjInRoi(imgVesselSkel, roi, tools.roiDilation);

                        tools.print("Computing distance between each cell and its closest vessel...");
                        tools.computeCellsVesselsDists(popCellRoi, imgVesselDistMapInv);
                    }
                    
                    tools.print("Drawing results...");
                    tools.drawResultsInRoi(popCellRoi, objVesselRoiDil, objSkelRoi, imhCell, imhCellDist, imhVessel, imhSkel);   
                    
                    tools.print("Writing results...");
                    tools.writeResultsInRoi(popCellRoi, objSkelRoiDil, objSkelRoi, imgVesselDistMap, imgCell, roi, vesselVolRoi, outDir, imgName);
                }
                
                tools.saveCloseDrawings(imgCell, imgVessel, imhCell, imhCellDist, imhVessel, imhSkel, outDir, imgName);
                tools.saveRois(rois, outDir, imgName);
                
                tools.closeImage(imgCell);
                tools.closeImage(imgCellMask);
                if (imgVessel != null){
                    tools.closeImage(imgVessel);
                    tools.closeImage(imgVesselMask);
                    tools.closeImage(imgVesselSkel);
                    tools.closeImage(imgVesselDistMap.getImagePlus());
                    tools.closeImage(imgVesselDistMapInv.getImagePlus());
                }
            }
            tools.closeResults();
            tools.print("--- All done! ---");
        } catch (DependencyException | ServiceException | IOException | FormatException ex) {
            Logger.getLogger(VeCell.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(VeCell.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}

