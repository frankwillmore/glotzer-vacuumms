/* pmft2pov.c */

//  IN:    x\ty\tz\tF12\n PMFT vals
//  OUT:   .pov  povray format

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <ftw_param.h>
#include <ftw_pov.h>


int main(int argc, char *argv[]) {

  char line[256];
  char *xs, *ys, *zs, *pmft_s;
  double x, y, z;
  char color[256];
  char phong_s[256];
  unsigned int pmft_color;
  unsigned int red, green, blue;
  double pmft;
  double transmittance = 0.7;
  double phong = 1.0;
  double blob_threshold = 0.65;
  double element_radius = 0.125;
  double box_x=10, box_y=10, box_z=10;
  double threshold_min = -1.0, threshold_max = 0.0; // range of values to display
  double pmft_color_min=-6.0, pmft_color_max=0.0; // range of values for color scale

  setCommandLineParameters(argc, argv);
  if (getFlagParam("--usage") || getFlagParam("-u")) {
    printf("usage:       pmft2pov    IN:      .pmf = potemntial of mean force (x, y, z, F12)\n");
    printf("                         OUT:     .pov = povray format (no header)\n\n");
    printf("\n");
    printf("                         --threshold_min [-1.0] \n");
    printf("                         --threshold_max [0.0] \n");
    printf("                         --pmft_color_min [-6.0] \n");
    printf("                         --pmft_color_max [0.0] \n");
    printf("                         --blob_threshold [0.65] \n");
    printf("                         --element_radius [0.125] \n");
    printf("                         --transmittance [0.7]  \n");
    printf("                         --phong [1.0] \n");
    printf("                         --blob \n");
    printf("                         --clip \n");
    printf("\n");
    exit(0);
  }
  writePOVHeader();

  getDoubleParam("--threshold_min", &threshold_min);
  getDoubleParam("--threshold_max", &threshold_max);
  getDoubleParam("--pmft_color_min", &pmft_color_min);
  getDoubleParam("--pmft_color_max", &pmft_color_max);
  getDoubleParam("--blob_threshold", &blob_threshold);
  getDoubleParam("--element_radius", &element_radius);
  getDoubleParam("--transmittance", &transmittance);
  getDoubleParam("--phong", &phong);

  printf("// begin pmft2pov records\n");
  
  if (getFlagParam("--clip")) printf("intersection { \n  box {<0,0,0>< %lf, %lf, %lf>} \n", box_x, box_y, box_z);
  if (getFlagParam("--blob")) printf("  blob { \n    threshold %lf \n", blob_threshold);

  // loop over all lines 
  while (1) {
    fgets(line, 256, stdin);
    if (feof(stdin)) break;
    
    xs = strtok(line, "\t");
    ys = strtok(NULL, "\t");
    zs = strtok(NULL, "\t");
    pmft_s = strtok(NULL, "\n");

    x = strtod(xs, NULL);
    y = strtod(ys, NULL);
    z = strtod(zs, NULL);
    pmft = strtod(pmft_s, NULL);
if (pmft < threshold_min || pmft > threshold_max) continue; // toss these

    // need to map pmft to a color, 0..511
    // pmft_color = floor(((pmft - pmft_color_min) / (pmft_color_max - pmft_color_min)) * 512);
    pmft_color = 511 - floor(((pmft - pmft_color_min) / (pmft_color_max - pmft_color_min)) * 512);

    if (pmft_color < 256) {
      red = 255 - pmft_color;
      green = pmft_color;
      blue = 0;
    }  
    else {
      red = 0;
      green = 511 - pmft_color;
      blue = pmft_color - 256;
    }

    sprintf (color, "rgb <%f, %f, %f>", red/256.0, green/256.0, blue/256.0);
    if (phong > 0) sprintf (phong_s, " finish {phong %lf} ", phong);
    else phong_s[0] = '0';

    //if (getFlagParam("--blob")) printf("    sphere{<%lf, %lf, %lf>, %lf %lf }\n", x, y, z, element_radius, 1.00000);
    if (getFlagParam("--blob")) printf("sphere{<%lf, %lf, %lf>, %lf 1.000000 texture{ pigment {color %s transmit %f} %s }}\n", x, y, z, element_radius, color, transmittance, phong_s);
    else {
      // We want transmittance to go to 1 as intensity goes to zero. 
      //transmittance = 1.0 - fv_intensity;
      transmittance = 0.0;
      printf("sphere{<%lf, %lf, %lf>, %lf texture{ pigment {color %s transmit %f} %s }}\n", x, y, z, element_radius, color, transmittance, phong_s);
    }
  } // end reading input

  if (getFlagParam("--blob") && getFlagParam("--clip")) {
     printf("  } // end blob \n");
     printf("  texture{\n");
     printf("    pigment {color %s transmit %lf } \n", color, transmittance);    
     printf("    finish {phong %lf} \n", phong);
     printf("  } // end texture \n");
     printf("} // end intersection \n");
  }
  else if (getFlagParam("--blob") && !getFlagParam("--clip")) {
//     printf("    texture{\n");
//     printf("      pigment {color %s transmit %lf } \n", color, transmittance);    
//     printf("      finish {phong %lf} \n", phong);
//     printf("    } // end texture \n");
     printf("  } // end blob \n");
  }
  else if (!getFlagParam("--blob") && getFlagParam("--clip")) {
     printf("} // end intersection \n");
  }
 
  fclose(stdin);

  printf("// end fvi2pov records\n");
  return 0;
}

