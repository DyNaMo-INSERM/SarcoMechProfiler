// directory to save CSV files and ROIs
saveDir = getDirectory("Choose a directory to save CSV and ROI files");
if (saveDir == "") {
    exit("No directory selected. Macro will exit.");
}

// Open the ROI Manager and clear previous ROIs
run("ROI Manager...");
roiManager("Reset");

// add multiple line ROIs interactively
do {
    waitForUser("Draw a line ROI and press OK to add it to the ROI Manager.");
    roiManager("Add"); // Add selected ROI to ROI Manager
    continueSelection = getBoolean("Would you like to add another ROI?");
} while (continueSelection == true);

// Get the number of ROIs in the ROI Manager
roiCount = roiManager("Count");
print("Total ROIs: " + roiCount); // Debugging log

if (roiCount == 0) {
    exit("No ROIs in ROI Manager. Macro will exit.");
}

// Get the list of ROI names
roiList = roiManager("List"); // This gives a string with all ROI names
roiNames = split(roiList, "\n"); // Split the string into an array based on newline

// Process each open image
openImages = getList("image.titles"); // Get the titles of all open images
for (imageIndex = 0; imageIndex < openImages.length; imageIndex++) {
    imageTitle = openImages[imageIndex];
    selectImage(imageTitle);  // Switch to the next image
    print("Processing image: " + imageTitle);

    // Choose the value column name based on the image title
    if (indexOf(toLowerCase(imageTitle), "beta") >= 0) {
        valueColumn = "beta";
    } else {
        valueColumn = "log10(E0)";
    }

    // Process each ROI for the current image
    for (i = 0; i < roiCount; i++) {
        roiManager("Select", i); // Select the ROI

        // Wait a little for the profile to be plotted and results to appear
        wait(500); // Give time for the results window to appear

        roiName = "ROI_" + i; // Generate a name for each ROI, can customize based on your needs
        print("Processing " + roiName); // Debugging log

        run("Clear Results"); 
        getPixelSize(unit, pixelWidth, pixelHeight);
        profile = getProfile();

        run("Clear Results"); 
        for (ip = 0; ip < profile.length; ip++) {
            setResult("distance_" + unit + "_", ip, ip * pixelHeight);
            setResult(valueColumn, ip, profile[ip]);
            updateResults();
        }

        // Create a unique file name for each image and ROI
        csvFilePath = saveDir + imageTitle + "_" + roiName + "_profile.csv";
        selectWindow("Results");

        // Save the profile as a CSV file
        saveAs("Results", csvFilePath);
    }
}

// Save all ROIs as a ZIP file
roiManager("Save", saveDir + "All_ROIs.zip");

// Show a completion message
showMessage("Processing Complete", "Profiles and ROIs saved to: " + saveDir);
