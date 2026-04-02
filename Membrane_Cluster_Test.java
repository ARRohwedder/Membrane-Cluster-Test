import ij.*;
import ij.plugin.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import java.util.*;
import ij.measure.*;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.*;
import ij.plugin.frame.PlugInFrame;
import ij.plugin.frame.RoiManager;


/** Dr. Arndt Rohwedder, Faculty of Biological Sciences, University of Leeds.
* The Membrane cluster plugin identifies the outer borders of cells and combined structures and analysis clusters on that border.
* It than analysis the internal particles in dependence to their position. Works for grey level images. */
public class Membrane_Cluster_Test implements PlugIn {

    static double newUpperLower = 0.2;
    static double newChannel1 = 1.0;
    static double newChannel2 = 2.0;
    static double newlarge = 6;
    static double newsmall = 0.1;

    public void run(String arg) {
        double pixres, membranelength, clustercount, clusterlength;
        int ch1, ch2;

        if (IJ.getImage().getBitDepth()!=8){IJ.showMessage("Gray scale measure", "Greyscale Image required");}

        ImagePlus imp1 = IJ.getImage();
        ImageStack stack = imp1.getStack();
        Calibration cal = imp1.getCalibration();
        ImageProcessor ip = (ImageProcessor)imp1.getProcessor();
        pixres = cal.pixelWidth;
        int size = stack.getSize();


//Dialog box
        GenericDialog gd = new GenericDialog("Threshold for detection:", IJ.getInstance());
        gd.addNumericField("Conf. interv. for co-localization: ", newUpperLower, 3);
        gd.addMessage("Select Channels to compare");
        gd.addNumericField("Channel 1", newChannel1,0);
        gd.addNumericField("Channel 2", newChannel2,0);
        gd.addMessage("Enter desired particle size (um2)");
        gd.addNumericField("Largest", newlarge,3);
        gd.addNumericField("Smallest", newsmall,3);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        if (gd.invalidNumber()) {
            IJ.showMessage("Error", "Threshold not in range");
            return;
        }
        newUpperLower = gd.getNextNumber();
        newChannel1 = gd.getNextNumber();
        ch1 = (int)newChannel1;
        newChannel2 = gd.getNextNumber();
        ch2 = (int)newChannel2;
        newlarge = gd.getNextNumber();
        double wurzel = Math.sqrt(newlarge);
        int large = (int) Math.round(Math.pow(wurzel/pixres,2));
        newsmall = gd.getNextNumber();
        wurzel = Math.sqrt(newsmall);
        int small = (int) Math.round(Math.pow(wurzel/pixres,2));



//measure image dimensions
        int w = ip.getWidth();
        int h = ip.getHeight();

        ImagePlus first = IJ.createImage("Outlines", "8-bit black", w, h, 1);
        ImagePlus clusterimp = IJ.createImage("Cluster", "8-bit black", w, h, 1);
        ImageProcessor firstip = (ImageProcessor)first.getProcessor();
        ImageProcessor clusterip = (ImageProcessor)clusterimp.getProcessor();
        ImageProcessor stackip = (ImageProcessor)imp1.getProcessor();
        //ImageProcessor stackworkip = stack.getProcessor(1);
        ImageProcessor stackworkip = stack.getProcessor(ch1);
        ImageProcessor workip = (ImageProcessor)stackworkip.duplicate();
        first.setCalibration(imp1.getCalibration());
        clusterimp.setCalibration(imp1.getCalibration());
        ResultsTable rt = new ResultsTable();
        ResultsTable rt3 = new ResultsTable();

//initialise variables with 0

        clustercount=0;
        membranelength = 0;
        clusterlength = 0;
        double greyvalue = 0;
//Detection procedure

        ImagePlus findcell = IJ.createImage("findcell", "8-bit black", w, h,1);
        ImageProcessor findcellip = (ImageProcessor)findcell.getProcessor();

//sum all slices and put them into dedicated image
        for (int width=1; width<=w-1; width++) {
                for (int height=1; height<=h-1; height++) {
                    for (int slice=1; slice<=size; slice++) {
                        stackip = stack.getProcessor(slice);
                        greyvalue = (stackip.get(width,height)*1.2)+greyvalue;
                    }
                    if (greyvalue > 254){greyvalue = 255;}
                    findcellip.putPixel(width,height,(int) Math.round(greyvalue));
                    greyvalue = 0;
                }
        }

//use threshold on cell image
        findcellip.invert();
        IJ.setAutoThreshold(findcell, "Huang");
        RankFilters desp = new RankFilters();
        desp.rank(findcellip, 1.0, desp.MEDIAN,desp.BRIGHT_OUTLIERS,0);

//identify cell outline and store coordinates, draw frame to avoid cell touching image border
        RoiManager managercell = RoiManager.getInstance();
        if (managercell == null){managercell = new RoiManager();}
        ResultsTable tempResultscell = new ResultsTable();
        ParticleAnalyzer partycell = new ParticleAnalyzer( ParticleAnalyzer.ADD_TO_MANAGER + ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES + ParticleAnalyzer.INCLUDE_HOLES,Measurements.CENTER_OF_MASS + Measurements.AREA,tempResultscell,((w*h/100)*25),((w-50)*(h-50)));
        Roi frame = new Roi(0,0,w,h);
        findcellip.setColor(0);
        frame.setStrokeWidth(2);
        findcellip.drawRoi(frame);
        partycell.analyze(findcell);

//select the largest identified cell
        int howmany = managercell.getCount();
        if (howmany<1){IJ.showMessage("Outline could not be identified");return;}
        int largest = 0;
        if (howmany > 1){
            for (int which = 1; which < howmany; which++){
                if (tempResultscell.getValue("Area",which)>largest){largest=which;}
            }
        }


//coordinates of centre of mass
        int roicenterx = (int) tempResultscell.getValue("XM",largest);
        int roicentery = (int) tempResultscell.getValue("YM",largest);

//transfer of identified cell to poly ROI and get coordinates for all points into array and paint it to image
        Roi cell = managercell.getRoi(largest);

        firstip.setColor(255);
        firstip.draw(cell);

//clean work copy for further identification
        workip.setColor(0);
        workip.fillOutside(cell);
        FloatPolygon cellfloat = cell.getFloatPolygon();
        int lang = cellfloat.npoints;
        int xcell [] = new int[lang];
        int ycell [] = new int[lang];
        for (int fill=0; fill < lang; fill++){
                xcell[fill] = (int) Math.round(cellfloat.xpoints[fill]);
                ycell[fill] = (int) Math.round(cellfloat.ypoints[fill]);
        }
        tempResultscell.reset();
        managercell.reset();


// distance (hypothenuse) of centre of mass and angle to centre is needed. Using circular calculation, stored in angle.
        double angle [] = new double[lang];
        double hypothenuse [] = new double[lang];
        double xdiff [] = new double [lang];
        double xdiffcast = 0;
        double xdiffsquare = 0;
        double ydiffsquare = 0;

//find angle and radius to centre of mass
        for (int position=0; position < lang; position++){
                xdiffsquare = Math.pow(Math.abs(roicentery-ycell[position]),2);
                ydiffsquare = Math.pow(Math.abs(roicenterx-xcell[position]),2);
                hypothenuse [position]= Math.sqrt(xdiffsquare+ydiffsquare);
                xdiffcast = cellfloat.xpoints[position];
                xdiff [position] = Math.abs(roicenterx-xdiffcast);
                angle[position] = Math.toDegrees(Math.asin(xdiff[position]/hypothenuse[position]));
        }

// calculate average grey values for selected membranes and delete pixels in workip.
        double membrane1 [] = new double[lang];
        double membrane2 [] = new double[lang];
        for (int position=0; position < lang; position++){
            stackip = stack.getProcessor(ch1);
            for (int pos=0; pos<4; pos++){
                membrane1[position] = membrane1[position] + stackip.get(((int) Math.round((Math.sin(Math.toRadians(angle[position]))*(hypothenuse[position]-pos))+roicenterx)),((int) Math.round((Math.cos(Math.toRadians(angle[position]))*(hypothenuse[position]-pos))+roicentery)));
                workip.putPixel(((int) Math.round((Math.sin(Math.toRadians(angle[position]))*(hypothenuse[position]-pos))+roicenterx)),((int) Math.round((Math.cos(Math.toRadians(angle[position]))*(hypothenuse[position]-pos))+roicentery)),0);
            }
            stackip = stack.getProcessor(ch2);
            for (int pos=0; pos<4; pos++){
                membrane2[position] = membrane2[position] + stackip.get(((int) Math.round((Math.sin(Math.toRadians(angle[position]))*(hypothenuse[position]-pos))+roicenterx)),((int) Math.round((Math.cos(Math.toRadians(angle[position]))*(hypothenuse[position]-pos))+roicentery)));
            }
            membrane1[position] = membrane1[position] / 4;
            membrane2[position] = membrane2[position] / 4;
        }

// normalize grey intensities for membranes in each channel
        double average1 = 0;
        double average2 = 0;
        for (int position=0; position < lang; position++){
            average1 = average1 + membrane1[position];
            average2 = average2 + membrane2[position];
        }
        average1 = average1/lang;
        if (average1<10){IJ.showMessage("Warning: average of reference intensity very low");}
        average2 = average2/lang;
        if (average2<10){IJ.showMessage("Warning: average of second image intensity very low");}

// find co-localization and store it in new array
        double membcolocal [] = new double [lang];
        for (int position=0; position < lang; position++){
            if ((membrane1[position]/average1)/(membrane2[position]/average2) > 1-newUpperLower && (membrane1[position]/average1)/(membrane2[position]/average2)< 1+newUpperLower){
                membcolocal [position] = (membrane1[position]/average1)/(membrane2[position]/average2);
            }
            else {membcolocal[position]=0;}

        }

// find clusters of co-localization and calculate size and distribution along membrane. Representation drawn to clusterip
        int cluster [] =  new int [lang];
        for (int position=0; position< lang-1; position++) {
                if (membcolocal[position]>0 && membcolocal[position+1]>0){
                    cluster[position] = 255;
                    if (Math.abs(membcolocal[position-1]-membcolocal[position])>0){
                            clustercount = clustercount+1;
                            }
                }
                while (membcolocal[position]>0 && membcolocal[position+1]>0){
                    clusterlength = Math.sqrt(Math.pow(Math.abs(xcell[position]-xcell[position+1]),2)+Math.pow(Math.abs(ycell[position]-ycell[position+1]),2))+clusterlength;
                    clusterip.putPixel(xcell[position],ycell[position],255);
                    position++;
                }
                if (clusterlength > 0){
                    rt3.incrementCounter();
                    rt3.addValue("cluster size", clusterlength*pixres);
                }
                clusterlength = 0;
                membranelength = Math.sqrt(Math.pow(Math.abs(xcell[position]-xcell[position+1]),2)+Math.pow(Math.abs(ycell[position]-ycell[position+1]),2))+membranelength;
        }

        clusterimp.show();
        rt3.show("Cluster sizes");
        IJ.log("Identified clusters: "+clustercount+"  Membrane length : "+(membranelength*pixres)+" Cluster per um "+(clustercount/(membranelength*pixres)));

//find vesicles and their co-localization in workip
        workip.invert();
        ResultsTable tempResults = new ResultsTable();
        RoiManager manager = RoiManager.getInstance();
        if (manager == null){
            manager = new RoiManager();
        }
//analyse internal particles in workip using threshold
        ParticleAnalyzer party = new ParticleAnalyzer( ParticleAnalyzer.ADD_TO_MANAGER,Measurements.CENTER_OF_MASS,tempResults,small,large);
        //stackip = stack.getProcessor(ch1);
        ImagePlus workimp = new ImagePlus("work",workip); //working imageplus copy

        IJ.setAutoThreshold(workimp, "Default"); //build binary
        party.analyze(workimp);
        int wieviel = manager.getCount();
        if (wieviel<1){IJ.showMessage("No vesicles found in outline");return;}
        double roix [] = new double[wieviel];
        double roiy [] = new double[wieviel];
//get positions of ROIs and draw them into firstip
        for (int roipos=1; roipos<wieviel-1;roipos++) {
                roix[roipos] = tempResults.getValue("XM",roipos);
                roiy[roipos] = tempResults.getValue("YM",roipos);
                firstip.setColor(155);
                firstip.draw(manager.getRoi(roipos));
        }

//measure grey values and size of ROIs
        tempResults = manager.multiMeasure(imp1);
        double roigreys1 [] = new double[wieviel];
        double roigreys2 [] = new double[wieviel];
        double roiareas [] = new double[wieviel];

        String wo;
        for (int roipos=1; roipos<wieviel-1;roipos++) {
                wo = "Area"+roipos;
            roiareas[roipos] = tempResults.getValue(wo, (ch1-1));
                wo = "Mean"+roipos;

            roigreys1[roipos] = tempResults.getValue(wo, (ch1-1));

            roigreys2[roipos] = tempResults.getValue(wo, (ch2-1));

        }

//calculate distance of identified vesicles to outer surface
        double membdist [] = new double[wieviel];
        double vector = 0;
        double roiaverage1 = 0;
        double roiaverage2 = 0;
        double vectorx;
        double vectory;

        for (int roipos=0; roipos<wieviel;roipos++) {
                membdist [roipos] = (double) w;
                for (int degree=0; degree<lang;degree++) {
                        vectorx = Math.abs(xcell[degree]-(roix[roipos]));
                        vectorx = Math.pow(vectorx,2);
                        vectory = Math.pow(Math.abs(ycell[degree]-(roiy[roipos])),2);
                        vector = Math.sqrt(vectorx+vectory);
                        if (vector<membdist[roipos]){membdist[roipos] = vector;}
                }
        }

//normalize grey values for each slice
        for (int roipos=1; roipos<wieviel;roipos++) {
                roiaverage1 = roiaverage1 + roigreys1[roipos];
                roiaverage2 = roiaverage2 + roigreys2[roipos];
        }
        roiaverage1 = roiaverage1/wieviel;
        roiaverage2 = roiaverage2/wieviel;
        double [] normgrey1 = new double[wieviel];
        double [] normgrey2 = new double[wieviel];
        for (int roipos=1; roipos<wieviel;roipos++) {
                normgrey1[roipos] = roigreys1[roipos]/roiaverage1;
                normgrey2[roipos] = roigreys2[roipos]/roiaverage2;
        }
//calculate co-localized vesicles and add all values to output table
        double [] colocal = new double[wieviel];
        for (int roipos=1; roipos<wieviel;roipos++) {
                if ((normgrey1[roipos]/normgrey2[roipos] > 1-newUpperLower) && (normgrey1[roipos]/normgrey2[roipos] < 1+newUpperLower)){colocal[roipos]=normgrey1[roipos]/normgrey2[roipos];}
                else {colocal[roipos]=0;}
        }
        ResultsTable roiResults = new ResultsTable();
        for (int roipos=1; roipos<wieviel;roipos++) {
            roiResults.incrementCounter();
            roiResults.addValue("x pos", roix[roipos]);
            roiResults.addValue("y pos", roiy[roipos]);
            roiResults.addValue("Area (um2)", roiareas[roipos]);
            roiResults.addValue("Membrane distance (um)", membdist[roipos]*pixres);
            roiResults.addValue("Colocalisation", colocal[roipos]);
        }
        roiResults.show("Vesicles");
        first.show();
        }
    }
