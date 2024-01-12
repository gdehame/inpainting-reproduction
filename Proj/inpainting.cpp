#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <math.h>
#include <unistd.h>
//Command-line parsing
#include <CLI11.hpp>

//Image filtering and I/O
#define cimg_display 0
#include <CImg.h>
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

// RENDU LATEX + CODE POUR LE 28/04 (faut absolument envoyer un mail, possibilité de juste dire "je rends la semaine pro" et le rendre la semaine suivante)


//Global flag to silent verbose messages
bool silent;
//Normalization factor for the data term
double alpha = 255;

int patchSize = 9;

using namespace std;

vector<int> gradientX(unsigned char * image, vector<int> & filter, int width, int height, int nbChannels)
{
  //Auxiliary function to calculate the horizontal gradients of the image
  //It converts the image in gray scale then compute gradients because it's only defined on gray scale
  vector<int> grad(width*height);
  for (int i = 0; i < width; i++)
  {
    for (int j = 0; j < height; j++)
    {
      int nextPixel = nbChannels*(width*j+((width + (i+1) % width) % width));
      int previousPixel = nbChannels*(width*j+((width + (i-1) % width) % width));
      int indPixel = nbChannels*(width*j+i);
      if (filter[width*j+((width + (i+1) % width) % width)] == 0)
      {
        //We do forward difference if possible
        grad[j*width+i] = ((int)image[nextPixel] - (int)image[indPixel] + (int)image[nextPixel+1] - (int)image[indPixel+1] + (int)image[nextPixel+2] - (int)image[indPixel+2] )/3;
      }
      else if (filter[width*j+((width + (i-1) % width) % width)] == 0)
      {
        //If not we do backward difference
        grad[j*width+i] = ((int)image[previousPixel] - (int)image[indPixel] + (int)image[previousPixel+1] - (int)image[indPixel+1] + (int)image[previousPixel+2] - (int)image[indPixel+2] )/3;
      }
      else
      {
        //If still not then we are not on the frontier of the mask so we don't care
        grad[j*width+i] = 0;
      }
    }
  }
  return grad;
}

vector<int> gradientY(unsigned char * image, vector<int>& filter, int width, int height, int nbChannels)
{
  //Auxiliary function to calculate the vertical gradients of the image
  //It converts the image in gray scale then compute gradients because it's only defined on gray scale
  vector<int> grad(width*height);
  for (int i = 0; i < width; i++)
  {
    for (int j = 0; j < height; j++)
    {
      int nextPixel = nbChannels*(width*((height + (j+1) % height) % height)+i);
      int previousPixel = nbChannels*(width*((height + (j-1) % height) % height)+i);
      int indPixel = nbChannels*(width*j+i);
      if (filter[width*((height + (j+1) % height) % height)+i] == 0)
      {
        //We do forward difference if possible
        grad[j*width+i] = ((int)image[nextPixel] - (int)image[indPixel] + (int)image[nextPixel+1] - (int)image[indPixel+1] + (int)image[nextPixel+2] - (int)image[indPixel+2] )/3;
      }
      else if (filter[width*((height + (j-1) % height) % height)+i] == 0)
      {
        //If not we do backward difference
        grad[j*width+i] = ((int)image[previousPixel] - (int)image[indPixel] + (int)image[previousPixel+1] - (int)image[indPixel+1] + (int)image[previousPixel+2] - (int)image[indPixel+2] )/3;
      }
      else
      {
        //If still not then we are not on the frontier of the mask so we don't care
        grad[j*width+i] = 0;
      }
    }
  }
  return grad;
}

void updateGradients(vector<unsigned char>& image, vector<int>& gradX, vector<int>& gradY, array<int,2>& center, vector<int> & filter, int width, int height, int nbChannels)
{
  //Auxiliary function to update the gradients
  //Litterally does the same as gradientX and gradientY but modifying the gradX and gradY parameters
  for (int i = center[0]-(patchSize-1)/2 -1; i <= center[0]+(patchSize-1)/2 +1; i++)
  {
    if (i < 0 || i >= width) continue;
    for (int j = center[1]-(patchSize-1)/2 -1; j <= center[1]+(patchSize-1)/2 +1; j++)
    {
      if (j < 0 || j >= height) continue;

      int nextPixel = nbChannels*(width*j+((width + (i+1) % width) % width));
      int previousPixel = nbChannels*(width*j+((width + (i-1) % width) % width));
      int indPixel = nbChannels*(width*j+i);
      if (filter[width*j+((width + (i+1) % width) % width)] == 0)
      {
        gradX[j*width+i] = ((int)image[nextPixel] - (int)image[indPixel] + (int)image[nextPixel+1] - (int)image[indPixel+1] + (int)image[nextPixel+2] - (int)image[indPixel+2] )/3;
      }
      else if (filter[width*j+((width + (i-1) % width) % width)] == 0)
      {
        gradX[j*width+i] = ((int)image[previousPixel] - (int)image[indPixel] + (int)image[previousPixel+1] - (int)image[indPixel+1] + (int)image[previousPixel+2] - (int)image[indPixel+2] )/3;
      }
      else
      {
        gradX[j*width+i] = 0;
      }

      nextPixel = nbChannels*(width*((height + (j+1) % height) % height)+i);
      previousPixel = nbChannels*(width*((height + (j-1) % height) % height)+i);
      if (filter[width*((height + (j+1) % height) % height)+i] == 0)
      {
        gradY[j*width+i] = ((int)image[nextPixel] - (int)image[indPixel] + (int)image[nextPixel+1] - (int)image[indPixel+1] + (int)image[nextPixel+2] - (int)image[indPixel+2] )/3;
      }
      else if (filter[width*((height + (j-1) % height) % height)+i] == 0)
      {
        gradY[j*width+i] = ((int)image[previousPixel] - (int)image[indPixel] + (int)image[previousPixel+1] - (int)image[indPixel+1] + (int)image[previousPixel+2] - (int)image[indPixel+2] )/3;
      }
      else
      {
        gradY[j*width+i] = 0;
      }
    }
  }
}

vector<array<double, 2>> normal(vector<int>& filter, int width, int height)
{
  //Compute normals of the mask
  vector<array<double,2>> normals(height*width);
  normals[0] = {-10,-10}; //used to check if it failed
  for (int i = 0; i < width; i++)
  {
    for (int j = 0; j < height; j++)
    {
      int upPixel = width*((height + (j+1) % height) % height)+i;
      int centerPixel = width*j+i;
      int rightPixel = width*j+((i+1) % width) % width;
      double nx = filter[rightPixel] - filter[centerPixel];
      double ny = filter[upPixel] - filter[centerPixel];
      double norm = sqrt(nx*nx+ny*ny);
      if (norm == 0)
      {
        normals[j*width+i] = {0,0};
      }
      else
      {
        normals[j*width+i] = {nx / norm, ny / norm};
      }
    }
  }
  return normals;
}

void updateNormals(vector<array<double, 2>>& normals, vector<int>& filter, array<int,2>& center, int width, int height)
{
  //Updates the normals, does the same as normal but modifying the parameter normals
  for (int i = center[0]-(patchSize-1)/2 -1; i <= center[0]+(patchSize-1)/2 +1; i++)
  {
    if (i < 0 || i >= width) continue;
    for (int j = center[1]-(patchSize-1)/2 -1; j <= center[1]+(patchSize-1)/2 +1; j++)
    {
      if (j < 0 || j >= height) continue;

      int upPixel = width*((height + (j+1) % height) % height)+i;
      int centerPixel = width*j+i;
      int rightPixel = width*j+((i+1) % width) % width;
      double nx = filter[rightPixel] - filter[centerPixel];
      double ny = filter[upPixel] - filter[centerPixel];
      double norm = sqrt(nx*nx+ny*ny);
      if (norm == 0)
      {
        normals[j*width+i] = {0,0};
      }
      else
      {
        normals[j*width+i] = {nx / norm, ny / norm};
      }
    }
  }
}

void updateConfidencesAndFilter(vector<double>& confidence, vector<int>& filter, array<int,2>& argMaxPriority, int width, int height)
{
  //Function to update confidence of the target patch we just replaced and remove them from the filter (they are not to be replaced anymore)
 	double newConfidence = 0;
	int cardPatch = 0;
  //First loops to compute the new confidence
  for (int i = -(patchSize-1)/2; i <= (patchSize-1)/2; i++)
  {
    if (i + argMaxPriority[0] < 0 || i + argMaxPriority[0] >= width)
    {
      continue;
    }
    for (int j = -(patchSize-1)/2; j <= (patchSize-1)/2; j++)
    {
      if (0 > j + argMaxPriority[1] || j + argMaxPriority[1] >= height)
      {
        continue;
      }
      cardPatch++;
      if (filter[width*(j + argMaxPriority[1]) + i + argMaxPriority[0]] == 0)
      {
        newConfidence += confidence[width*(j + argMaxPriority[1]) + i + argMaxPriority[0]];
      }      
    }
  }
  newConfidence /= (double) cardPatch;
  //Second ones to update the values
  for (int i = -(patchSize-1)/2; i <= (patchSize-1)/2; i++)
  {
    if (i + argMaxPriority[0] < 0 || i + argMaxPriority[0] >= width)
    {
      continue;
    }
    for (int j = -(patchSize-1)/2; j <= (patchSize-1)/2; j++)
    {
      if (0 > j + argMaxPriority[1] || j + argMaxPriority[1] >= height)
      {
        continue;
      }
      if (filter[width*(j + argMaxPriority[1]) + i + argMaxPriority[0]] == 1)
      {
        confidence[width*(j + argMaxPriority[1]) + i + argMaxPriority[0]] = newConfidence;
        filter[width*(j + argMaxPriority[1]) + i + argMaxPriority[0]] = 0;
      }
    }
  } 
}

void copyPatch(vector<unsigned char>& output, vector<int>& filter, array<int,2>& argMaxPriority, array<int,2>& closestPatchCenter, int width, int height, int nbChannels)
{
  //Copies the patch whose center is closestPatchCenter in the one whose center is argMaxPriority in output
  for (int i = -(patchSize-1)/2; i <= (patchSize-1)/2; i++)
  {
    if (i + argMaxPriority[0] < 0 || i + argMaxPriority[0] >= width)
    {
      continue;
    }
    for (int j = -(patchSize-1)/2; j <= (patchSize-1)/2; j++)
    {
      if (0 > j + argMaxPriority[1] || j + argMaxPriority[1] >= height || filter[width * (j + argMaxPriority[1]) + i + argMaxPriority[0]] == 0)
      {
        continue;
      }
      int indexPixel = nbChannels*(width * (j + argMaxPriority[1]) + i + argMaxPriority[0]);
      int indexSrcPixel = nbChannels*(width * (j + closestPatchCenter[1]) + i + closestPatchCenter[0]);
      output[indexPixel] = output[indexSrcPixel];
      output[indexPixel+1] = output[indexSrcPixel+1];
      output[indexPixel+2] = output[indexSrcPixel+2];
    }
  }
}

void updateFrontier(vector<array<int,2>>& frontier, vector<array<double,2>>& normals,  array<int,2>& argMaxPriority, int width, int height)
{
  //Updates the frontier using the normals
  int xUpperBound = argMaxPriority[0] + (patchSize-1)/2;
  int xLowerBound = argMaxPriority[0] - (patchSize-1)/2;
  int yUpperBound = argMaxPriority[1] + (patchSize-1)/2;
  int yLowerBound = argMaxPriority[1] - (patchSize-1)/2;
  frontier.erase(remove_if(frontier.begin(), frontier.end(), [xLowerBound, xUpperBound, yUpperBound, yLowerBound](array<int,2> p) {return (p[0] >= xLowerBound && p[0] <= xUpperBound && p[1] >= yLowerBound && p[1] <= yUpperBound);}), frontier.end());
  for (int i = xLowerBound -1; i <= xUpperBound +1; i++)
  {
    if (i < 0 || i >= width) continue;
    for (int j = yLowerBound -1; j <= yUpperBound +1; j++)
    {
      if (j < 0 || j >= height) continue;
      if (normals[j*width+i][0] != 0 || normals[j*width+i][1] != 0)
      { //Then it's part of the frontier
        frontier.push_back({i,j});
      }
    }
  }
}

bool patchInSrcZone(vector<int>& filter, int x, int y, int width, int height)
{
  //Checks if the patch centered in (x,y) is entirely in the source zone or not
  for (int i = -(patchSize-1)/2; i <= (patchSize-1)/2; i++)
  {
    if (x+i < 0 || x+i >= width) continue;
    for (int j = -(patchSize-1)/2; j <= (patchSize-1)/2; j++)
    {
      if (y+j < 0 || y+j >= height) continue;
      if (filter[width*(y+j)+x+i] == 1)
      {
        return false;
      }
    }
  }
  return true;
}

array<int,3> convertRGBToCIELab(unsigned char r, unsigned char g, unsigned char b)
{
  //Translate RGB to CIE Lab based on a conversion I found on internet
  //I have no idea of what is CIE Lab nor the logic of that conversion
  array<int,3> rgbPixel = {(int)r,(int)g,(int)b};
  float fr = (float)r;
  float fg = (float)g;
  float fb = (float)b;
  array<int,3> cieLabPixel = {0,0,0};
  float x = (0.41253*fr + 0.357580*fg + 0.180423*fb)/0.950456;
  float y = 0.212671*fr + 0.71516*fg + 0.072169*fb;
  float z = (0.019334*fr + 0.119193*fg + 0.950227*fb)/1.088754;
  float l = 0;
  float fy = 0;
  if (y > 0.008856)
  {
    fy = cbrt(y);
    l = 116.0*fy - 16.0;
  }
  else
  {
    fy = 7.787*x+16.0/116.0;
    l = 903.3*y;
  }
  float fx = 0;
  float fz = 0;
  if (x > 0.008856)
  {
    fx = cbrt(x);
  }
  else
  {
    fx = 7.787*x+16.0/116.0;
  }
  if (z > 0.008856)
  {
    fz = cbrt(z);
  }
  else
  {
    fz = 7.787*z+16.0/116.0;
  }
  int li = (int) (l*255.0/100.0);
  int a = (int) 500.0*(fx-fy);
  int br = (int) 200.0*(fy-fz);
  //cout << "blem: a=" << a << " et b=" << br << "(" << 500.0*(fx-fy) << "," << 200.0*(fy-fz) << ") et x=" << x << " et y=" << y << " et z=" << z << " et fx=" << fx << " et fy=" << fy << " et fz=" << fz << "   " << (int)r << "," << (int)g << "," << (int)b << endl;
  cieLabPixel[0] = li;
  cieLabPixel[1] = a;
  cieLabPixel[2] = br;
  return cieLabPixel;
}

int main(int argc, char **argv)
{
  //First part is just the usual loader you gave us in TP
  CLI::App app{"inpainting"};
  std::string sourceImage;
  app.add_option("-s,--source", sourceImage, "Source image")->required()->check(CLI::ExistingFile);;
  std::string targetImage;
  app.add_option("-m,--mask", targetImage, "Mask image")->required()->check(CLI::ExistingFile);;
  std::string outputImage= "output.png";
  app.add_option("-o,--output", outputImage, "Output image")->required();
  std::string outputMask= "./Proj/outputMask.png";
  bool stepByStepOutput = false;
  app.add_flag("--sbs", stepByStepOutput, "Step by step");
  silent = false;
  app.add_flag("--silent", silent, "No verbose messages");
  CLI11_PARSE(app, argc, argv);

  //Image loading
  int width,height, nbChannels;
  unsigned char *source = stbi_load(sourceImage.c_str(), &width, &height, &nbChannels, 0);
  if (!silent) std::cout<< "Source image: "<<width<<"x"<<height<<"   ("<<nbChannels<<")"<< std::endl;
  int width_target,height_target, nbChannels_target;
  unsigned char *target = stbi_load(targetImage.c_str(), &width_target, &height_target, &nbChannels_target, 0);
  if (!silent) std::cout<< "Target image: "<<width_target<<"x"<<height_target<<"   ("<<nbChannels_target<<")"<< std::endl;
  if ((width*height) != (width_target*height_target))
  {
    std::cout<< "Image sizes do not match. "<<std::endl;
    exit(1);
  }
  if (nbChannels < 3)
  {
    std::cout<< "Input images must be RGB images."<<std::endl;
    exit(1);
  }
  if (nbChannels_target > 4)
  {
    cout << "Mask must be a black and white or RGB image" << endl; 
    exit(1);
  }


  //Main computation
  std::vector<unsigned char> output(width*height*nbChannels);
  for(auto i = 0 ; i < width ; ++i)
  {
    for(auto j = 0; j < height; ++j)
    {
      auto indexPixel = nbChannels*(width*j+i);
      unsigned char r = source[ indexPixel ];
      unsigned char g = source[ indexPixel + 1];
      unsigned char b = source[ indexPixel + 2];
      output[ indexPixel ] = r;
      output[ indexPixel + 1 ] = g;
      output[ indexPixel + 2 ] = b;
      if (nbChannels == 4) //just copying the alpha value if any
        output[ indexPixel + 3] = source[ indexPixel + 3];
    }
  }

  vector<int> filter(width*height);
  //We initialize the filter
  for (int i = 0; i < width; i++)
  {
    for (int j = 0; j < height; j++)
    {
      int indexPixel = nbChannels_target*(width*j+i);
      if (nbChannels_target >= 3)
      {
        if ((int) target[indexPixel] <= 50  && (int) target[indexPixel+1] <= 50 && (int) target[indexPixel+2] <= 50)
        {
          filter[j*width+i] = 0;
        }
        else
        {
          filter[j*width+i] = 1;
        }
      }
      else
      {
        if ((int) target[indexPixel] == 0)
        {
          filter[j*width+i] = 0;
        }
        else if ((int) target[indexPixel] == 255)
        {
          filter[j*width+i] = 1;
        }
        else
        {
          cout << "Error: target image should be black and white (it's the filter) " << (int)target[indexPixel] << "," << (int) target[indexPixel+1] << endl;
          exit(1);
        }
      }
    }
  }

  vector<int> gradX = gradientX(source, filter, width, height, nbChannels);
  vector<int> gradY = gradientY(source, filter, width, height, nbChannels);  

  //Normals are the gradients of the target filter, thus they are 0 outside the frontier
  vector<array<double,2>> normals = normal(filter, width, height);

  vector<double> confidence(width*height);
  //We initialize the confidence
  for (int i = 0; i < width; i++)
  {
    for (int j = 0; j < height; j++)
    {
      confidence[j*width+i] = 1 - filter[j*width+i];
    }
  }

  vector<array<int, 2>> frontier;
  //Building the frontier
  for (int i = 0; i < width; i++)
  {
    for (int j = 0; j < height; j++)
    {
      if (normals[j*width+i][0] != 0 || normals[j*width+i][1] != 0)
      { //Then it's part of the frontier
        frontier.push_back({i,j});
      }
    }
  }

  int nbIter = 0;
  while (frontier.size() > 0)
  {
    //Calculating priorities and priority maximizer
    vector<double> priorities(frontier.size());
    double maxPriority = -1;
    array<int,2> argMaxPriority = {-1,-1};
    for (int i = 0; i < priorities.size(); i++)
    {
      int pixelX = frontier[i][0];
      int pixelY = frontier[i][1];
      double conf = 0;
      int card = 0;
      for (int j = max(0, pixelX - (patchSize-1)/2); j <= min(width-1, pixelX + (patchSize-1)/2); j++)
      {
        for(int k = max(0, pixelY - (patchSize-1)/2); k <= min(height-1, pixelY + (patchSize-1)/2); k++)
        {
          card++;
          conf += confidence[width*k+j];
        }
      }
      conf /= (double) card;
      priorities[i] = conf * abs(-gradY[pixelY*width+pixelX] * normals[pixelY*width+pixelX][0] + gradX[pixelY*width+pixelX] * normals[pixelY*width+pixelX][1]) / alpha;
      if (i == 0 || maxPriority < priorities[i])
      {
        maxPriority = priorities[i];
        argMaxPriority = frontier[i];
      }
    }

    //Looking for closestPatch
    array<int, 2> closestPatchCenter;
    int minDistance = -1;
    for (int i = min((patchSize-1)/2,argMaxPriority[0]); i <= max(width - (patchSize-1)/2 -1, argMaxPriority[0]); i++)
    {
    	for (int j = min((patchSize-1)/2, argMaxPriority[1]); j <= max(height - (patchSize-1)/2 -1, argMaxPriority[1]); j++)
      {
        if (!patchInSrcZone(filter, i, j, width, height))
        {
          //We only copy patches that are entirely in the source region
          continue;
        }
        int distance = 0;
        int distEucl = sqrt((i - argMaxPriority[0])*(i - argMaxPriority[0])+(j - argMaxPriority[1])*(j - argMaxPriority[1]));
        bool calculated = false; //Used to distinguish between distance 0 because exactly the same and because no pixel were compared
        //We iterate over all the patch
        for (int k = -(patchSize-1)/2; k <= (patchSize-1)/2; k++)
        {
          if (minDistance != -1 && distance >= minDistance) break;

          if (argMaxPriority[0] + k < 0 || argMaxPriority[0] + k >= width)
          {
            //We ignore points out of the image
            continue;
          }
          for (int l = -(patchSize-1)/2; l <= (patchSize-1)/2; l++)
          {
            if (minDistance != -1 && distance >= minDistance) break;

            if (argMaxPriority[1] + l < 0 || argMaxPriority[1] + l >= height || filter[width*(argMaxPriority[1]+l)+argMaxPriority[0]+k] == 1)
            {
              //We ignore points out of the image
              continue;
            }
            calculated = true;
            int indexPatchPixel = nbChannels*(width*(argMaxPriority[1] + l) + argMaxPriority[0] + k);
						int indexSrcPixel = nbChannels*(width*(j+l)+i+k);
            array<int,3> patchPixel = convertRGBToCIELab(output[indexPatchPixel], output[indexPatchPixel+1], output[indexPatchPixel+2]);
            array<int,3> srcPixel = convertRGBToCIELab(output[indexSrcPixel], output[indexSrcPixel+1], output[indexSrcPixel+2]);
            for (int m = 0; m < 3; m++)
            {
              distance += (srcPixel[m] - patchPixel[m])*(srcPixel[m] - patchPixel[m]);
            }
          }
        }
				if (minDistance == -1 || (calculated && distance + distEucl < minDistance))
				{
          //We add a distance between coordinates of the centers of patches in the image to break equalities
					minDistance = distance + distEucl;
					closestPatchCenter = {i,j};
				}
      }
    }

    //Copying data from source patch to most prioritized patch
    copyPatch(output, filter, argMaxPriority, closestPatchCenter, width, height, nbChannels);
    
    //Updating confidence, filter, normals and gradients then frontier for the next loop iteration
		updateConfidencesAndFilter(confidence, filter, argMaxPriority, width, height);
    updateNormals(normals, filter, argMaxPriority, width, height);
    updateGradients(output, gradX, gradY, argMaxPriority, filter, width, height, nbChannels);
    updateFrontier(frontier, normals, argMaxPriority, width, height);

    if (stepByStepOutput)
    {
      //For debugging purposes, outputs the image at each step
      stbi_write_png(outputImage.c_str(), width, height, nbChannels, output.data(), nbChannels*width);

      sleep(1);

      vector<unsigned char> sbsMask(3*width*height);
      for (int i = 0; i < width; i++)
      {
        for (int j = 0; j < height; j++)
        {
          int indexPxl = nbChannels*(j*width+i);
          sbsMask[indexPxl] = 255*filter[j*width+i];
          sbsMask[indexPxl+1] = 255*filter[j*width+i];
          sbsMask[indexPxl+2] = 255*filter[j*width+i];
        }
      }
      stbi_write_png(outputMask.c_str(), width, height, 3, sbsMask.data(), 3*width);
    }
  }


  //Final export
  if (!silent) std::cout<<"Exporting.."<<std::endl;
  int errcode = stbi_write_png(outputImage.c_str(), width, height, nbChannels, output.data(), nbChannels*width);
  if (!errcode)
  {
    std::cout<<"Error while exporting the resulting image. Error code : " << errcode <<std::endl;
    exit(errcode);
  }

  stbi_image_free(source);
  stbi_image_free(target);
  exit(0);
}
