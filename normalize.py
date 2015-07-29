# SCRIPT WRITTEN BY JOSEPH ANTONAKAKIS
# Reorienting atoms above a sheet based upon atomic coordinates and 
# sheet given in original coordinate system.

# INPUT FILE - 
# First 5 lines: 
# Coordinates of sheet in 1st coordinate system. 
# Next 5 lines (NO SPACE IN BETWEEN SETS): 
# Coordinates of sheet in 2nd coordinate system. 
# Lines 11 on (NO SPACE IN BETWEEN SETS): 
# Coordinates of atoms relative to 1st sheet in 1st coordinate system.
# Final space: 
# After all lines of data to indicate no more data. 
# NOTE: 
# FIRST SHEET'S TWO ATOMS' COORDINATES MUST BE ANALAGOUS TO THE 
# SECOND SHEET'S TWO ATOMS' COORDINATES!



from sympy import *
import math

### HELPER METHODS ###

# Accepts: 3 floats 
# Returns: the square root of the sum of the squares of 3 numbers 
# USED AS DIVIDING FACTOR FOR VECTOR NORMALIZATION 
def squareRootOf3(a, b, c): 
  return math.sqrt(a**2 + b**2 + c**2)


# Builds X and Z matrices for 3D linear regression
# Accepts: Empty arrays and modifies them accordingly
# Xarr results in an array with the following pattern:
# [1 x1 y1]
# [1 x2 y2]
# [1 x3 y3]
# ...
# Zarr results in an array with the following pattern: 
# [z1]
# [z2]
# [z3]
# ... 
# MEANT TO BE USED ITERATIVELY
def buildXandZArrays(Xarr, Zarr, columns):
  innerArray = [1]
  mat = Matrix() 
  for x in range(0, len(columns)-1):
    innerArray.append(float(columns[x]))
  Xarr.append(innerArray)
  Zarr.append([columns[len(columns)-1]])


# Accepts: X and Z matrices built for 3D linear regression
# Returns: BMatrix with coefficients for a best-fit plane eqaution
# for a give set of coordinates in 3D of the following form: 
# [B0]
# [B1]
# [B2] --> B0 + B1*x + B2*y = z
def findBMat(XMat, ZMat):
  BMat = (XMat.transpose()*XMat).inv()*XMat.transpose()*ZMat 
  return BMat


# Accepts: BMatrix that was the result of a 3D linear regression
# Returns: String form of the plane equation 
def printPlaneEq(BMat): 
  return ("" + str(float(BMat[0])) + " + " + str(float(BMat[1])) + "*x + " 
  + str(float(BMat[2])) + "*y = z")


# Accepts: BMatrix that was the result of a 3D linear regression
# Returns: Matrix object that is the plane's normal vector
def findNormal(BMatrix): 
  dividingFactor = squareRootOf3(BMatrix[1], BMatrix[2], -1)
  array = [ float(BMatrix[1])/dividingFactor, 
            float(BMatrix[2])/dividingFactor, 
            -1/dividingFactor ]
  return Matrix(array)


# Accepts: 2 vectors as matrices in 3 dimensions
# Returns: Matrix object that is the cross product of those parameters
def planeNormalCrossProduct(vector1, vector2): 
  resultArray = []
  resultArray.append((vector1[1]*vector2[2])-(vector1[2]*vector2[1]))
  resultArray.append(-(vector1[0]*vector2[2])+(vector1[2]*vector2[0]))
  resultArray.append((vector1[0]*vector2[1])-(vector1[1]*vector2[0]))
  return Matrix(resultArray)


# Accepts: 2 orthogonal vectors and uses matrix row reduction to 
# determine a 3rd orthogonal vector 
# Returns: Matrix object that represents that 3rd orthogonal vector
def findThirdOrthogVector(vector1, vector2): 
  mat = vector1.transpose().col_join(vector2.transpose()).rref()
  if mat[1][1] == 1: 
    dividingFactor = squareRootOf3((mat[0][2]), (mat[0][5]), (1))
    resultArray = [ -mat[0][2]/dividingFactor, 
                    -mat[0][5]/dividingFactor, 
                    1/dividingFactor ]

    return Matrix(resultArray)
  elif mat[1][1] == 2:
    dividingFactor = squareRootOf3((mat[0][1]), 1, (mat[0][4]))
    resultArray =[ -mat[0][1]/dividingFactor, 
                   1/dividingFactor,
                   -mat[0][4]/dividingFactor ]
    return Matrix(resultArray)


# Accepts: X and Z matrices from two directions of a plane in a coord system
# Returns: Array that represents a basis for that coord system 
def createCoordSystem(X5, Z5): 
  B5 = findBMat(X5, Z5)
  planeNormal5 = findNormal(B5)
  dX = X5[4] - X5[1] # 2nd point x minus 1st point x
  dY = X5[5] - X5[2] # 2nd point y minus 1st point y 
  z2ndPoint = B5[0] + B5[1]*X5[4] + B5[2]*X5[5]
  z1stPoint = B5[0] + B5[1]*X5[1] + B5[2]*X5[2]
  dZ = z2ndPoint - z1stPoint
  dividingFactor = squareRootOf3(dX, dY, dZ)  
  planeOrthog = Matrix([[dX/dividingFactor], 
                        [dY/dividingFactor], 
                        [dZ/dividingFactor]])
  # planeOrthog is a vector that's always in a KNOWN orientation 
  thirdOrthog = findThirdOrthogVector(planeOrthog, planeNormal5)
  return [planeOrthog, thirdOrthog, planeNormal5]


# Accepts: Array of 3 vectors (Matrix objects) that are in 3 dimensions 
# Returns: Returns 3x3 Matrix that combines said 3 vectors 
def complete3By3Matrix(vectorArray):
  assert isinstance(vectorArray, list)
  assert len(vectorArray) == 3
  return vectorArray[0].row_join(vectorArray[1]).row_join(vectorArray[2])


# Accepts: Two arrays that consist of Matrix objects that represent a basis 
# of a coordinate system
# Returns: Change of basis matrix calcuated in terms of the second basis array
# (1st basis --> 2nd basis)
# Not used to project the molecules, but rather express in terms of the second coord sys
def changeOfBasisMatrix(basis1, basis2):
  assert isinstance(basis1, list)
  assert isinstance(basis2, list)
  basis2Matrix = complete3By3Matrix(basis2)
  resultingVectors = []
  for x in basis1: 
    resultingVectors.append(basis2Matrix.inv()*x)
  finalMat = complete3By3Matrix(resultingVectors)
  return finalMat


# Accepts: An array that indicates a 3D coordinate system basis 
# Prints out the cross products of 1st vector dot 2nd, 2nd dot 3rd, 
# and 1st dot 3rd 
# USED FOR TESTING PURPOSES 
# Iteration could have been used but I was lazy 
def printDotProducts(firstCoordSys): 
  oneDotTwo = (firstCoordSys[0][0]*firstCoordSys[1][0] + 
    firstCoordSys[0][1]*firstCoordSys[1][1] + 
    firstCoordSys[0][2]*firstCoordSys[1][2]) 
  twoDotThree = (firstCoordSys[1][0]*firstCoordSys[2][0] + 
    firstCoordSys[1][1]*firstCoordSys[2][1] + 
    firstCoordSys[1][2]*firstCoordSys[2][2])
  oneDotThree = (firstCoordSys[2][0]*firstCoordSys[0][0] + 
    firstCoordSys[2][1]*firstCoordSys[0][1] + 
    firstCoordSys[2][2]*firstCoordSys[0][2])

  print oneDotTwo
  print twoDotThree 
  print oneDotThree

# Accepts: Array of arrays that have 3 coordinates per element 
# Returns: A list of the points outputed properly for testing with 
# the following website: http://www.bodurov.com/VectorVisualizer
def printPoint(element): 
  outputted_string = "("
  for x in range(0, len(element)): 
    outputted_string += str(element[x])
    if x != len(element)-1:
      outputted_string += ","
  outputted_string += ":"
  for x in range(0, len(element)): 
    outputted_string += str(element[x])
    if x != len(element)-1:
      outputted_string += ","
  outputted_string += ")"
  print outputted_string


def printVector(element): 
  outputted_string = "("
  for x in range(0, len(element)):
    outputted_string += str(element[x])
    if x != len(element)-1:
      outputted_string += ","
  outputted_string += ")"
  print outputted_string


  # Used for the purposes of displaying the elements' coordinates
  def dispFinalAtoms(vectorArray, coordSys, offsetMatrix, a, b, origin):
  assert a == 1 or a == -1
  assert b == 1 or b == -1
  for x in vectorArray:
    assert len(x) == 3
    finalVector = (x[2][0]*coordSys[0] + x[2][1]*a*coordSys[1] +
    + x[2][2]*b*coordSys[2])
    finalVector += (origin - 
                    offsetMatrix[0]*coordSys[0] -
                    offsetMatrix[1]*a*coordSys[1] - 
                    offsetMatrix[2]*b*coordSys[2])
    print (x[0] + "   " + str(finalVector[0]) + "   " + 
      str(finalVector[1]) + "   " + str(finalVector[2]))


### PART ONE: FIND TWO PLANES PER STRUCTURE ###

# First sheet's 5 point regression
X1_5 = Matrix(); X1_5Arr = []; Z1_5 = Matrix(); Z1_5Arr = []
# Second sheet's 5 point regression
X2_5 = Matrix(); X2_5Arr = []; Z2_5 = Matrix(); Z2_5Arr = []

# List of other vectors
vectors = []

allSurfacePoints = []


loopVar = 0
f = open("sampleData.txt")
loopVar = 0
for line in f:
  columns = line.split()
  if len(columns) == 0: # Reached the end
    break
  # First 5 
  if loopVar < 5:
    innerArray = []
    for x in range(0, len(columns)):
      innerArray.append(float(columns[x]))
    allSurfacePoints.append(innerArray)
    buildXandZArrays(X1_5Arr, Z1_5Arr, columns)
  # Next 5
  elif loopVar >= 5 and loopVar < 10:
    innerArray = []
    for x in range(0, len(columns)):
        innerArray.append(float(columns[x]))
    allSurfacePoints.append(innerArray)
    buildXandZArrays(X2_5Arr, Z2_5Arr, columns)
  # Creating a list of vectors 
  elif loopVar >= 10: 
    assert len(columns) == 4 # Symbol, x, y, z
    newVector = Matrix([[columns[1]], [columns[2]], [columns[3]]])
    vectors.append([columns[0], newVector]) # Column's first column is symbol 
    # NOTE: Creating an array b/c going to add coefficients later
  # Iterating through loop variable 
  loopVar += 1
f.close() 



X1_5 = Matrix(X1_5Arr)
Z1_5 = Matrix(Z1_5Arr)

X2_5 = Matrix(X2_5Arr)
Z2_5 = Matrix(Z2_5Arr)

firstB5 = findBMat(X1_5, Z1_5)
secondB5 = findBMat(X2_5, Z2_5)

# FLAG: FIRST COORDINATE SYSTEM BASIS ARRAY
firstCoordSys = createCoordSystem(X1_5, Z1_5)



# FLAG: SECOND COORDINATE SYSTEM BASIS ARRAY
secondCoordSys = createCoordSystem(X2_5, Z2_5)
secondCoordSys[1] = -secondCoordSys[1]



### PART TWO: FIND CHANGE OF BASIS MATRIX ###
# Finding the change of basis matrix 
changeMat = changeOfBasisMatrix(firstCoordSys, secondCoordSys)


### PART THREE: FIND EACH ATOM IN TERMS OF 1ST COORD SYSTEM (COEFFICIENTS) ###

firstCoordMatrix = complete3By3Matrix(firstCoordSys)
for x in vectors: 
  resultingVector = firstCoordMatrix.inv()*x[1]
  x.append(resultingVector) 
  # In terms of 1st coordinate system basis vectors 




### PART FOUR: PRINT OUT EACH VALUE AND COORD'S IN TERMS OF 2ND COORD SYSTEM 
### RELATIVE TO SECOND SHEET ###

origin = allSurfacePoints[0]
offsetCoefficients = firstCoordMatrix.inv()*Matrix(origin)
# Offset coefficients * corresponding basis vector must be added 
# to each coordinate of the sorbed atoms 
origin2 = allSurfacePoints[5]
origin2 = Matrix(origin2)



print "FINAL ATOMS - FIRST CONFIG"
# Essentially, using the coefficients of the point RELATIVE TO 1ST SHEET
# and multiplying them by the analagous basis vectors of the 2ND SHEET
dispFinalAtoms(vectors, secondCoordSys, offsetCoefficients, 1, 1, origin2)
print "FINAL ATOMS - SECOND CONFIG"
dispFinalAtoms(vectors, secondCoordSys, offsetCoefficients, 1, -1, origin2)
print "FINAL ATOMS - THIRD CONFIG"
dispFinalAtoms(vectors, secondCoordSys, offsetCoefficients, -1, 1, origin2)
print "FINAL ATOMS - FOURTH CONFIG"
dispFinalAtoms(vectors, secondCoordSys, offsetCoefficients, -1, -1, origin2)







