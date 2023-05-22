# Imports
import cv2
import numpy as np

def areaFilter(minArea, inputImage):
    # Perform an area filter on the binary blobs:
    componentsNumber, labeledImage, componentStats, componentCentroids = \
        cv2.connectedComponentsWithStats(inputImage, connectivity=4)

    # Get the indices/labels of the remaining components based on the area stat
    # (skip the background component at index 0)
    remainingComponentLabels = [i for i in range(1, componentsNumber) if componentStats[i][4] >= minArea]

    # Filter the labeled pixels based on the remaining labels,
    # assign pixel intensity to 255 (uint8) for the remaining pixels
    filteredImage = np.where(np.isin(labeledImage, remainingComponentLabels) == True, 255, 0).astype('uint8')

    return filteredImage

# Read image
imageFile = "vanRees2011DNS_vorticity.png"
inputImage = cv2.imread(imageFile)

# Conversion to CMYK (just the K channel):

# Convert to float and divide by 255:
imgFloat = inputImage.astype(np.float) / 255.

# Calculate channel K:
kChannel = 1 - np.max(imgFloat, axis=2)

# Convert back to uint 8:
kChannel = (255 * kChannel).astype(np.uint8)

# Threshold image:
binaryThresh = 190
_, binaryImage = cv2.threshold(kChannel, binaryThresh, 255, cv2.THRESH_BINARY)


# Filter small blobs:
minArea = 100
binaryImage = areaFilter(minArea, binaryImage)

# Use a little bit of morphology to clean the mask:
# Set kernel (structuring element) size:
kernelSize = 3
# Set morph operation iterations:
opIterations = 2
# Get the structuring element:
morphKernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernelSize, kernelSize))
# Perform closing:
binaryImage = cv2.morphologyEx(binaryImage, cv2.MORPH_CLOSE, morphKernel, None, None, opIterations, cv2.BORDER_REFLECT101)
cv2.imwrite("vanRees2011DNS_vorticity-only_black.png",binaryImage)

# cv2.imshow("binaryImage [closed]", binaryImage)
# cv2.waitKey(0)

