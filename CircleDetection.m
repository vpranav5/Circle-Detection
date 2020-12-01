% Name: Pranav Teja Varanasi
% UT EID: ptv247

% hough transform

% detect edge using standard method

twins = imread('twins.jpg');

snake = imread('snake.jpg');

planets = imread('planets.jpg');

gumballs = imread('gumballs.jpg');

coins = imread('coins.jpg');

sr71 = imread('sr71.jpg');

solar = imread('solar.jpg');

% figure;
% imshow(twins);
% 
% figure;
% imshow(snake);
% 
% figure;
% imshow(planets);
% 
% figure;
% imshow(gumballs);
% 
% figure;
% imshow(coins);



houghCenters = detectCirclesHT(gumballs, 60);

ransacCenters = detectCirclesRANSAC(planets, 100);

clusterMatrix = clusterPixels(sr71, 10);

% returns an n x 2 matrix with the coordinates of all the points considered
% to be boundaries

boundaryMatrix = boundaryPixels(clusterMatrix);

markedImage = markImage(sr71, boundaryMatrix);

function houghCenters = detectCirclesHT(image, radius) % Input is color image
   
   binPercent = 0.0001;
    
   grayScale = rgb2gray(image);
   smoothImage = imgaussfilt(grayScale);
   
   % Use canny edge detector on smoothed image for greater accuracy
   detectEdges = edge(smoothImage, 'Canny'); 
    
   numRows = size(detectEdges, 1);
   numCols = size(detectEdges, 2);
   
   % 2 pi radian constant
   twoPI = 2 * pi;
   
   houghMatrix = zeros(numRows, numCols);
   
   for row = 1:numRows % loop through every edge pixel (x, y)
       for col = 1:numCols
           % loop through all the points of a circle
           edgeCenter = detectEdges(row, col);
            
           if(edgeCenter == 1)
               % draw a circle
               % make a result matrix which is copy of image to store frequencies
               % go thru all the points on the circle and increment one
            
               for radianAngle = 1:0.001:twoPI
                    estimateA = round(row + radius * cos(radianAngle));
                    estimateB = round(col - radius * sin(radianAngle));
                    
                    if(estimateA >= 1 && estimateA <= numRows)
                        if(estimateB >= 1 && estimateB <= numCols)
                            houghMatrix(estimateA, estimateB) = houghMatrix(estimateA, estimateB) + 1;  
                        end
                        
                    end
                    
               end
               
               
           end  
           
       end 
   end
   
   %figure;
   imshow(image);
      
   
    houghRows = size(houghMatrix, 1);
    houghCols = size(houghMatrix, 2);
    
    totalElements = houghRows * houghCols;
    topTen = round(binPercent * totalElements);
     
    [sortedMatrix originalIndices] = sort(houghMatrix(:), 'descend'); %convert the matrix to a vector and sort it
    [houghRow houghCol] = ind2sub(size(houghMatrix), originalIndices(1:topTen));
     
    % houghRow and houghCol are matrices
    rowBound = size(houghRow, 1);
    
    centersMatrix = zeros(rowBound, 2); % n x 2 matrix of centers
    
     for index=1:rowBound
        checkCenter = [houghCol(index, 1), houghRow(index, 1)];
        viscircles(checkCenter, radius, 'Color', 'b');
        
        centersMatrix(index, 1) = houghRow(index, 1);
        centersMatrix(index, 2) = houghCol(index, 1);
        
     end
     
     % return the centers combined into a 2 x 1 matrix
     
    %imagesc(houghMatrix);
   
    houghCenters = centersMatrix;
    
    
end

function ransacCenters = detectCirclesRANSAC(image, radius) % Input is color image
    iterations = 1;
    binPercent = 0.1;

    thresholdDistance = 2;

    grayScale = rgb2gray(image);
    smoothImage = imgaussfilt(grayScale);
   
    % Use canny edge detector on smoothed image for greater accuracy
    detectEdges = edge(smoothImage, 'Canny'); 
    
    numRows = size(image, 1);
    numCols = size(image, 2);
    
    edgeRow = [];
    edgeCol = [];
    
    detectEdgeRows = size(detectEdges, 1);
    detectEdgeCols = size(detectEdges, 2);
    
    
    for edgeRowIndex = 1:detectEdgeRows
        
        for edgeColIndex = 1:detectEdgeCols
            
            edgeValue = detectEdges(edgeRowIndex, edgeColIndex);
            
            if(edgeValue == 1)
                
                edgeRow = [edgeRow; edgeRowIndex];
                edgeCol = [edgeCol; edgeColIndex];
                
            end
           
        end
        
    end
    
    voteMatrix = zeros(numRows, numCols);

    numEdgeRows = size(edgeRow, 1);
    
    numEdgeCols = size(edgeCol, 1);
    
    
     for iteration = 1:iterations


         % pick 3 random edge points, if the edge points form a line diagonally, vertically, horizontally, 
         % or if the calculated radius is greater then input radius then keep looking for points
         % setup system of equations using bisector method and pass it to matlab syseq function

         % Bisector Method Source: https://www.qc.edu.hk/math/Advanced%20Level/circle%20given%203%20points.htm


         invalidPoints = true;

         centerRow = 0;
         centerCol = 0;
         centerRadius = 0;

         while(invalidPoints)

            randRowOne = randi(numEdgeRows);
            randColOne = randi(numEdgeCols); 

            randRowTwo = randi(numEdgeRows);
            randColTwo = randi(numEdgeCols); 

            randRowThree = randi(numEdgeRows);
            randColThree = randi(numEdgeCols); 


            firstPointRow = edgeRow(randRowOne);
            firstPointCol = edgeCol(randColOne);

            secondPointRow = edgeRow(randRowTwo);
            secondPointCol = edgeCol(randColTwo);

            thirdPointRow = edgeRow(randRowThree);
            thirdPointCol = edgeCol(randColThree);


            % if the 3 points together form a line diagonally, vertically, or horizontally, skip to next iteratioon


            % points form a horizontal line
            if(firstPointRow == secondPointRow && secondPointRow == thirdPointRow)
                "skipping horizontal";
                continue;
            end

            % points form a vertical line since row represents y and col represents x
            if(firstPointCol == secondPointCol && secondPointCol == thirdPointCol)
                "skipping vertical";
                
                % go to next loop iteration and dont try to fit points
                continue;
            end

            slopeOneTwo =  (secondPointRow - firstPointRow) / (secondPointCol - firstPointCol) ;

            slopeTwoThree = (thirdPointRow - secondPointRow) / (thirdPointCol - secondPointCol) ;

            if(slopeOneTwo == slopeTwoThree)
                "skipping diagonal";
                continue;
            end

            "reaching circle calculation";
            % if this code executes we have found three valid points

            % setup system of equations using bisector method to calculate center of circle
            % perpendicular bisectors of two chords meet at the center of circle


            midOneTwoRow = (firstPointRow + secondPointRow) / 2;
            midOneTwoCol = (firstPointCol + secondPointCol) / 2;

            midTwoThreeRow = (secondPointRow + thirdPointRow) / 2;
            midTwoThreeCol = (secondPointCol + thirdPointCol) / 2;

            % bisector L1 is perpendicular to PQ, take negative inverse

            bisectorOneSlope = -1 / slopeOneTwo;

            % bisector L2 is perpendicular to QR

            bisectorTwoSlope = -1 / slopeTwoThree;

            % setup system of equations for bisectors using form (y - y1) = m (x - x1)

            % declare a system of equations in terms of x and y
            syms x y;

            % row represents the "y" coordinate and col represents "x" coordinate
            
            % move variables to left side and constants to equation right side

            % setup equation of first bisector
            % y - midOneTwoRow == bisectorOneSlope * (x - midOneTwoCol) ;

            equationOne = (-1 * bisectorOneSlope * x) + y == (midOneTwoRow + (-1 * bisectorOneSlope * midOneTwoCol));

            % setup equation of second bisector
            % y2 - midTwoThreeRow == bisectorTwoSlope * (x2 - midTwoThreeCol) ;

            equationTwo = (-1 * bisectorTwoSlope * x) + y == (midTwoThreeRow + (-1 * bisectorTwoSlope * midTwoThreeCol));

            % solve system of equations setup as AX = B

            [xCoff,yCoff] = equationsToMatrix([equationOne, equationTwo], [x, y]);
            
            
            % check if there is a solution to system of equations
            
            % count the number of NaN or Inf values in result matrices
            numberXNan = sum(isnan(xCoff(:)));
            numberXInf = sum(isinf(xCoff(:)));
            
            numberYNan = sum(isnan(yCoff(:)));
            numberYInf = sum(isinf(yCoff(:)));
            
            % skip to next iteration if no solution to system of equations
            if(numberXNan > 0 || numberXInf > 0 || numberYNan > 0 || numberYInf > 0)
                continue;
            end

            centerSolution = linsolve(xCoff,yCoff);


            % x value represents a column value based on inputs to equation
            % cast to an unsigned int to get a valid image coordinate
            
            calculateCol = uint32(centerSolution(1)); 

            % y value represents a row value based on inputs to equation
            calculateRow = uint32(centerSolution(2));
            
            % don't consider circles with center outside image
            if(calculateCol < 1 || calculateCol > numCols || calculateRow < 1 || calculateRow > numRows)
                continue;
            end


            % the radius for this new circle is the distance between a point, arbitrarily pick point 3, and the new center

            doubleCalcCol = double(calculateCol);

            diffCol = doubleCalcCol - thirdPointCol;

            doubleCalcRow = double(calculateRow);

            diffRow = doubleCalcRow - thirdPointRow;

           
            squareDist = (diffCol) ^ 2 + (diffRow) ^ 2;

            squareDist = double(squareDist);

            calcRadius = sqrt(squareDist);

            if(calcRadius <= radius)
                centerCol = calculateCol;
                centerRow = calculateRow;
                centerRadius = calcRadius;
                "reaching here";
                
                invalidPoints = false; % stop searching for points as we have found a candidate circle

            % end radius if statement
            end


         % end while loop
         end
     
     % end outer for loop
     end

    % row and col of the candidate circle
    % these are valid row and col values in original image
    centerRow;
    centerCol;
    
    
     % go through all the edge points and find all the edge points that are within a certain threshold distance to and from this center point

     % if the distance is within a threshold outside and inside th boundary of circle, add the number of inlier pixels found for this circle center to vote voteMatrix
    
     % form "doughnut" area around circle for where inliers can exist
     % inliers can be found in radius + - thresholdDistance
     
     outerBound = centerRadius + thresholdDistance;
     
     innerBound = centerRadius - thresholdDistance;
     
     for edgeIndex = 1:numEdgeRows
         
         currentEdgeRow = double(edgeRow(edgeIndex));
         currentEdgeCol = double(edgeCol(edgeIndex));
         
         % get distance between predicted center and the current edge point
         rowDiff = currentEdgeRow - centerRow;
         colDiff = currentEdgeCol - centerCol;
         
         squareEdgeDist = (rowDiff) ^ 2 + (colDiff) ^ 2;
         
         squareEdgeDist = double(squareEdgeDist);
         
         edgeCenterDist = sqrt(squareEdgeDist);

         % this edge's distance is valid for an inlier
         if(edgeCenterDist >= innerBound && edgeCenterDist <= outerBound)
             " reachine here";
             voteMatrix(centerRow, centerCol) = voteMatrix(centerRow, centerCol) + 1;
         end
         
     end
     
    voteRows = size(voteMatrix, 1);
    voteCols = size(voteMatrix, 2);
    
    totalElements = voteRows * voteCols;
    topPercent = round(binPercent * totalElements);
     
    [sortedMatrix originalIndices] = sort(voteMatrix(:), 'descend'); %convert the matrix to a vector and sort it
    
    [centerRow centerCol] = ind2sub(size(voteMatrix), originalIndices(1:topPercent));


    % centerRow and centerCol are matrices
    rowBound = size(centerRow, 1);
    
    centersMatrix = zeros(rowBound, 2); % n x 2 matrix of centers
    
    figure;
    imshow(image);
    
    for index=1:rowBound
        checkCenter = [centerCol(index, 1), centerRow(index, 1)];
        viscircles(checkCenter, radius, 'Color', 'b');
        
        centersMatrix(index, 1) = centerRow(index, 1);
        centersMatrix(index, 2) = centerCol(index, 1);
        
    end
    
    voteMatrix(centerRow(1), centerCol(1))
     
     
    % return the centers combined into a 2 x 1 matrix
     
    % imagesc(houghMatrix);
     
    ransacCenters = centersMatrix;
    
    
end



function clusterMatrix = clusterPixels(imageMatrix, kCenters)
    % imageMatrix

   numRows = size(imageMatrix, 1);
   numCols = size(imageMatrix, 2);
   
   % store which center each pixel belongs to, each pixel has a value associated between 1 to kCenters
   resultMatrix = zeros(numRows, numCols);


   % k x 4 matrix where 1st column is r value of center, 2nd column is g value of center, 3rd column is b value of center, 4th column is number of pixels in this center
   % need to clear 4th column in modifyCenters on each run
   % 4th column in originalCenters should always be 0
   
   % one matrix to keep track of previous center, and one matrix to keep track of the new centers, check if the matrices are equivalent
  

   originalCenters = zeros(kCenters, 4);
   modifyCenters = zeros(kCenters, 4);

   for index = 1:kCenters
       randRow = randi(numRows);
       randCol = randi(numCols);
       
       randRed = imageMatrix(randRow, randCol, 1);
       randGreen = imageMatrix(randRow, randCol, 2);
       randBlue = imageMatrix(randRow, randCol, 3);

       originalCenters(index, 1) = randRed;
       originalCenters(index, 2) = randGreen;
       originalCenters(index, 3) = randBlue;
       
       % 1 coordinate has been accumulated 
       originalCenters(index, 4) = 1;
       
   end
   
   modifyCenters = originalCenters;  % matrices are hard copies
  
   
   modifyCenters;
    
   centerChange = true; 
   
   % continue looping until centers converge to a value
   
   while(centerChange)
    
   
    % loop thru each pixel in image 

    for row = 1:numRows
       
        for col = 1:numCols
            
            pixelRed = double(imageMatrix(row, col, 1));
            pixelGreen = double(imageMatrix(row, col, 2));
            pixelBlue = double(imageMatrix(row, col, 3));
            
            
            minDistance = 1000000; 
           
            % track which center this pixel is closest to in terms of RGB values
            closestCenter = 0;
            
            
            for center = 1:kCenters
                
                 centerRed = originalCenters(center, 1);
                 centerGreen = originalCenters(center, 2);
                 centerBlue = originalCenters(center, 3);

                 
                 redDiff = (pixelRed - centerRed);
                 greenDiff = (pixelGreen - centerGreen);
                 blueDiff = (pixelBlue - centerBlue);
                
                 squareDistance = (redDiff ^ 2) + (greenDiff ^ 2) + (blueDiff ^ 2);
                 squareDistance = cast(squareDistance, 'double');
  
                 currentDistance = sqrt(squareDistance);
                 
                 if(currentDistance < minDistance)
                   
                     minDistance = currentDistance;
                  
                     closestCenter = center;
                 end
                 
            end

            % update result matrix with which center this pixel is closest to
            resultMatrix(row, col) = closestCenter;
               
            modifyCenters(closestCenter, 4) = modifyCenters(closestCenter, 4) + 1; % increment the 4th column storing how many pixels are associated with this center

            % add this pixels rgb values to the existing centers rgb values so we can take average later
            % cast to double before adding rgb values together

            modifyCenters(closestCenter, 1) = modifyCenters(closestCenter, 1) + pixelRed;
            modifyCenters(closestCenter, 2) = modifyCenters(closestCenter, 2) + pixelGreen;
            modifyCenters(closestCenter, 3) = modifyCenters(closestCenter, 3) + pixelBlue;
           
           
        % end inner for loop
        end
        
    % end outer for loop
    end
    
    modifyCenters;

    % take average of all pixels in modifyCenters to calculate new center

    for centerIndex = 1:kCenters

    	numPixels = modifyCenters(centerIndex, 4);
        
        if(numPixels > 0) % dont take average for centers that have no pixels in them
        
            modifyCenters(centerIndex, 1) = round(modifyCenters(centerIndex, 1) / numPixels);

            modifyCenters(centerIndex, 2) = round(modifyCenters(centerIndex, 2) / numPixels);

            modifyCenters(centerIndex, 3) = round(modifyCenters(centerIndex, 3) / numPixels);

        end 
        
    end
    
    modifyCenters;

    % set all the elements in the 4th column in modifyCenters to 1 since
    % only 1 element accumulated

    for loopCenters = 1:kCenters
    	modifyCenters(loopCenters, 4) = 1;
    end
    
    originalCenters;
    
    modifyCenters;

    % check if the originalCenters and modifyCenters matrices are equal 
    centerSame = isequal(originalCenters, modifyCenters);


    % if matrices are equal, then centerChange = false

    if (centerSame == true)
                
    	centerChange = false;
    else
    	% matrices are not equal, continue iterations to find a better center
        
    	% set the original centers matrix to modify centers 
        
    	originalCenters = modifyCenters;
    end
    
  % end outer while loop
  end


 % return result matrix as output 
 
 clusterMatrix = resultMatrix;

end


function boundaryMatrix = boundaryPixels(clusterMatrix)

       % image with boundaries marked
       % if the pixel next to one you are currently at is a different
       % cluster change the color
       
       % check up down left right
       
       numRows = size(clusterMatrix, 1);
       numCols = size(clusterMatrix, 2);
        
       boundaryRows = [];
       boundaryCols = [];
       
       for row = 1:numRows
           
           for col = 1:numCols


               % store which cluster from 1 to k this pixel belongs to 

               pixelCluster = clusterMatrix(row, col);
               
               
               % row and col change up, down, diag left, right moves
               rowChange = [-1,1,0,0,-1,-1,1,1];
               colChange = [0,0,-1,1,-1, 1,-1,1];
               
               numDifferent = 0;
               
               % check 4 directions from this pixel to see if it is on a boundary

               for index = 1:8
                   
                   moveRow = row + rowChange(index);
                   moveCol = col + colChange(index);
                   
                   if(moveRow >= 1 && moveRow <= numRows && moveCol >= 1 && moveCol <= numCols)
                       
                       moveCluster = clusterMatrix(moveRow, moveCol);
                       
                       if(moveCluster ~= pixelCluster)
                           numDifferent = numDifferent + 1;
                       end
                       
                       
                   end
                  
               end
               
               numDifferent;
               
               if(numDifferent > 0)
               	   % if this pixel is on a boundary, add it to the vectors keeping track of coordinates

                   boundaryRows = [boundaryRows; row];
               	   boundaryCols = [boundaryCols; col];
               end
              

           % end inner for loop
           end
       
       % end outer for loop  
       end
       
       resultMatrix = [boundaryRows boundaryCols];
        

       boundaryMatrix = resultMatrix;
       
       
end

function markedImage = markImage(imageMatrix, boundaryMatrix)
   % method to go thru the image matrix and mark all the pixels that are boundaries with yellow

   boundaryMatrix;
   
   % get number of entries vertically in boundary vectors

   numBoundRows = size(boundaryMatrix, 1);
   

   % loop thru the boundary matrix and change the color of coordinates in image matrix accordingly.

   for boundIndex = 1:numBoundRows

   		boundRow = boundaryMatrix(boundIndex, 1);
   		boundCol = boundaryMatrix(boundIndex, 2);

   		% color pixel yellow to show the boundary
   		imageMatrix(boundRow, boundCol, :) = [255,255,0];


   end
  
  figure;
  imshow(imageMatrix);
  markedImage = imageMatrix;

end 

 
 
 
 