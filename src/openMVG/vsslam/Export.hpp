#pragma once

#include <openMVG/third_party/png/png.h>
/*
    save_png
    -----------------------------------------------------

    A simple to save a png with a bit more flexibility. This function
    returns 0 on success otherwise -1.

    - filename:   the path where you want to save the png.
    - width:      width of the image
    - height:     height of the image
    - bitdepth:   how many bits per pixel (e.g. 8).
    - colortype:  PNG_COLOR_TYEP_GRAY
                  PNG_COLOR_TYPE_PALETTE
                  PNG_COLOR_TYPE_RGB
                  PNG_COLOR_TYPE_RGB_ALPHA
                  PNG_COLOR_TYPE_GRAY_ALPHA
                  PNG_COLOR_TYPE_RGBA          (alias for _RGB_ALPHA)
                  PNG_COLOR_TYPE_GA            (alias for _GRAY_ALPHA)
    - pitch:      The stride (e.g. '4 * width' for RGBA).
    - transform:  PNG_TRANSFORM_IDENTITY
                  PNG_TRANSFORM_PACKING
                  PNG_TRANSFORM_PACKSWAP
                  PNG_TRANSFORM_INVERT_MONO
                  PNG_TRANSFORM_SHIFT
                  PNG_TRANSFORM_BGR
                  PNG_TRANSFORM_SWAP_ALPHA
                  PNG_TRANSFORM_SWAP_ENDIAN
                  PNG_TRANSFORM_INVERT_ALPHA
                  PNG_TRANSFORM_STRIP_FILLER

 */
 static int save_png(std::string filename,
                     int width, int height, int bitdepth, int colortype,
                     unsigned char* data, int pitch,
                     int transform = PNG_TRANSFORM_IDENTITY);



 /* ----------------------------------------------------------------------- */

 static int save_png(std::string filename, int width, int height,
                     int bitdepth, int colortype,
                     unsigned char* data, int pitch, int transform)
 {
   int i = 0;
   int r = 0;
   FILE* fp = NULL;
   png_structp png_ptr = NULL;
   png_infop info_ptr = NULL;
   png_bytep* row_pointers = NULL;

   if (NULL == data) {
     printf("Error: failed to save the png because the given data is NULL.\n");
     r = -1;
     goto error;
   }

   if (0 == filename.size()) {
     printf("Error: failed to save the png because the given filename length is 0.\n");
     r = -2;
     goto error;
   }

   if (0 == pitch) {
     printf("Error: failed to save the png because the given pitch is 0.\n");
     r = -3;
     goto error;
   }

   fp = fopen(filename.c_str(), "wb");
   if (NULL == fp) {
     printf("Error: failed to open the png file: %s\n", filename.c_str());
     r = -4;
     goto error;
   }

   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
   if (NULL == png_ptr) {
     printf("Error: failed to create the png write struct.\n");
     r = -5;
     goto error;
   }

   info_ptr = png_create_info_struct(png_ptr);
   if (NULL == info_ptr) {
     printf("Error: failed to create the png info struct.\n");
     r = -6;
     goto error;
   }

   png_set_IHDR(png_ptr,
                info_ptr,
                width,
                height,
                bitdepth,                 /* e.g. 8 */
                colortype,                /* PNG_COLOR_TYPE_{GRAY, PALETTE, RGB, RGB_ALPHA, GRAY_ALPHA, RGBA, GA} */
                PNG_INTERLACE_NONE,       /* PNG_INTERLACE_{NONE, ADAM7 } */
                PNG_COMPRESSION_TYPE_BASE,
                PNG_FILTER_TYPE_BASE);

   row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);

   for (i = 0; i < height; ++i) {
     row_pointers[i] = data + i * pitch;
   }

   png_init_io(png_ptr, fp);
   png_set_rows(png_ptr, info_ptr, row_pointers);
   png_write_png(png_ptr, info_ptr, transform, NULL);

  error:

   if (NULL != fp) {
     fclose(fp);
     fp = NULL;
   }

   if (NULL != png_ptr) {

     if (NULL == info_ptr) {
       printf("Error: info ptr is null. not supposed to happen here.\n");
     }

     png_destroy_write_struct(&png_ptr, &info_ptr);
     png_ptr = NULL;
     info_ptr = NULL;
   }

   if (NULL != row_pointers) {
     free(row_pointers);
     row_pointers = NULL;
   }

   printf("And we're all free.\n");

   return r;
 }
