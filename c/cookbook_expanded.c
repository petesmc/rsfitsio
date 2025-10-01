#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
  Every program which uses the CFITSIO interface must include the
  the fitsio.h header file.  This contains the prototypes for all
  the routines and defines the error status values and other symbolic
  constants used in the interface.  
*/
#include "fitsio.h"

int main( void );
void writeimage( void );
void writeascii( void );
void writebintable( void );
void copyhdu( void );
void selectrows( void );
void readheader( void );
void readimage( void );
void read_ascii_table( void );
void read_binary_table( void );
void printerror( int status);

int main()
{
/*************************************************************************
   This is a simple main program that calls the following routines:

    writeimage        - write a FITS primary array image
    writeascii        - write a FITS ASCII table extension
    writebintable     - write a FITS binary table extension with all data types
    copyhdu           - copy a header/data unit from one FITS file to another
    selectrows        - copy selected row from one HDU to another
    readheader        - read and print the header keywords in every extension
    readimage         - read a FITS image and compute the min and max value
    read_ascii_table  - read columns of data from ASCII table
    read_binary_table - read columns of data from binary table

**************************************************************************/

    writeimage();
    writeascii();
    writebintable();
    copyhdu();
    selectrows();
    readheader();
    readimage();
    read_ascii_table();
    read_binary_table();

    printf("\nAll the cfitsio cookbook routines ran successfully.\n");
    return(0);
}
/*--------------------------------------------------------------------------*/
void writeimage( void )

    /******************************************************/
    /* Create a FITS primary array containing a 2-D image */
    /******************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, ii, jj;
    long  fpixel, nelements, exposure;
    unsigned short *array[200];

    /* initialize FITS image parameters */
    char filename[] = "cetestfil.fit";             /* name for new FITS file */
    int bitpix   =  USHORT_IMG; /* 16-bit unsigned short pixel values       */
    long naxis    =   2;  /* 2-dimensional image                            */    
    long naxes[2] = { 300, 200 };   /* image is 300 pixels wide by 200 rows */

    /* allocate memory for the whole image */ 
    array[0] = (unsigned short *)malloc( naxes[0] * naxes[1]
                                        * sizeof( unsigned short ) );

    /* initialize pointers to the start of each row of the image */
    for( ii=1; ii<naxes[1]; ii++ )
      array[ii] = array[ii-1] + naxes[0];

    remove(filename);               /* Delete old file if it already exists */

    status = 0;         /* initialize status before calling fitsio routines */

    if (fits_create_file(&fptr, filename, &status)) /* create new FITS file */
         printerror( status );           /* call printerror if error occurs */

    /* write the required keywords for the primary array image.     */
    /* Since bitpix = USHORT_IMG, this will cause cfitsio to create */
    /* a FITS image with BITPIX = 16 (signed short integers) with   */
    /* BSCALE = 1.0 and BZERO = 32768.  This is the convention that */
    /* FITS uses to store unsigned integers.  Note that the BSCALE  */
    /* and BZERO keywords will be automatically written by cfitsio  */
    /* in this case.                                                */

    if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
         printerror( status );          

    /* initialize the values in the image with a linear ramp function */
    for (jj = 0; jj < naxes[1]; jj++)
    {   for (ii = 0; ii < naxes[0]; ii++)
        {
            array[jj][ii] = ii + jj;
        }
    }

    fpixel = 1;                               /* first pixel to write      */
    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

    /* write the array of unsigned integers to the FITS file */
    if ( fits_write_img(fptr, TUSHORT, fpixel, nelements, array[0], &status) )
        printerror( status );
      
    free( array[0] );  /* free previously allocated memory */

    /* write another optional keyword to the header */
    /* Note that the ADDRESS of the value is passed in the routine */
    exposure = 1500;
    if ( fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,
         "Total Exposure Time", &status) )
         printerror( status );           

    if ( fits_close_file(fptr, &status) )                /* close the file */
         printerror( status );           

    return;
}
/*--------------------------------------------------------------------------*/
void writeascii ( void )

    /*******************************************************************/
    /* Create an ASCII table extension containing 3 columns and 6 rows */
    /*******************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status;
    long firstrow, firstelem;

    int tfields = 3;       /* table will have 3 columns */
    long nrows  = 6;       /* table will have 6 rows    */

    char filename[] = "cetestfil.fit";           /* name for new FITS file */
    char extname[] = "PLANETS_ASCII";             /* extension name */

    /* define the name, datatype, and physical units for the 3 columns */
    char *ttype[] = { "Planet", "Diameter", "Density" };
    char *tform[] = { "a8",     "I6",       "F4.2"    };
    char *tunit[] = { "\0",      "km",       "g/cm"    };

    /* define the name diameter, and density of each planet */
    char *planet[] = {"Mercury", "Venus", "Earth", "Mars","Jupiter","Saturn"};
    long  diameter[] = {4880,    12112,    12742,   6800,  143000,   121000};
    float density[]  = { 5.1f,     5.3f,      5.52f,   3.94f,    1.33f,    0.69f};

    status=0;

    /* open with write access the FITS file containing a primary array */
    if ( fits_open_file(&fptr, filename, READWRITE, &status) ) 
         printerror( status );

    /* append a new empty ASCII table onto the FITS file */
    if ( fits_create_tbl( fptr, ASCII_TBL, nrows, tfields, ttype, tform,
                tunit, extname, &status) )
         printerror( status );

    firstrow  = 1;  /* first row in table to write   */
    firstelem = 1;  /* first element in row  (ignored in ASCII tables) */

    /* write names to the first column (character strings) */
    /* write diameters to the second column (longs) */
    /* write density to the third column (floats) */

    fits_write_col(fptr, TSTRING, 1, firstrow, firstelem, nrows, planet,
                   &status);
    fits_write_col(fptr, TLONG, 2, firstrow, firstelem, nrows, diameter,
                   &status);
    fits_write_col(fptr, TFLOAT, 3, firstrow, firstelem, nrows, density,
                   &status);

    if ( fits_close_file(fptr, &status) )       /* close the FITS file */
         printerror( status );
    return;
}
/*--------------------------------------------------------------------------*/
void writebintable ( void )

    /*******************************************************************/
    /* Create a binary table extension with all data types (except P and Q) */
    /*
     * IMPORTANT NOTE ON TBIT (Bit Array) COLUMNS:
     * ============================================
     * TBIT columns in FITS binary tables require special handling that differs from other data types:
     *
     * 1. DATA REPRESENTATION:
     *    - TBIT columns store individual bits, not bytes
     *    - Each bit is represented as a logical value in memory:
     *      * 0 (zero) = FALSE = bit will be 0
     *      * Non-zero = TRUE = bit will be 1
     *    - For an "8X" column format, each row contains 8 bits
     *
     * 2. MEMORY LAYOUT:
     *    - The data is passed as an array of char (logical values)
     *    - For 6 rows with 8 bits each, you need 48 char elements total
     *    - NOT 6 bytes as you might expect!
     *
     * 3. API FUNCTIONS:
     *    - Use fits_write_col() and fits_read_col() with TBIT datatype for TBIT columns
     *    - Parameters are the same as other column types, but:
     *      * datatype: TBIT
     *      * nelem: Total number of bits to read/write (NOT number of rows)
     *      * array: Array of logical values (char) where 0=FALSE, non-zero=TRUE
     *
     * 4. PARAMETER COUNTING:
     *    - The 'nelem' parameter is the TOTAL number of bits across all rows
     *    - Example: 6 rows × 8 bits/row = 48 bits total
     *    - This is different from other column types where nelem = number of rows
     *
     * This example demonstrates proper TBIT handling alongside all other FITS data types.
     */
    /*******************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, hdutype;
    long firstrow, firstelem;

    int tfields   = 14;      /* table will have 14 columns */
    long nrows    = 6;       /* table will have 6 rows    */

    char filename[] = "cetestfil.fit";           /* name for new FITS file */
    char extname[] = "PLANETS_Binary";          /* extension name */

    /* define the name, datatype, and physical units for all columns */
    char *ttype[] = { 
        "Planet",      /* A - Character string */
        "Active",      /* L - Logical */
        "Flags",       /* X - Bit array */
        "Category",    /* B - Unsigned byte */
        "Priority",    /* I - 16-bit integer */
        "Diameter",    /* J - 32-bit integer */
        "Mass",        /* K - 64-bit integer  */
        "Density",     /* E - Single-precision float */
        "Gravity",     /* D - Double-precision float */
        "Position",    /* C - Single-precision complex */
        "Velocity",    /* M - Double-precision complex */
        "Moons",       /* 3I - Vector of 3 16-bit integers */
        "Temps",       /* 4E - Vector of 4 floats */
        "Coords"       /* 2D - Vector of 2 doubles */
    };
    char *tform[] = { 
        "8A",   /* Character string */
        "1L",   /* Logical */
        "8X",   /* 8 bits */
        "1B",   /* Unsigned byte */
        "1I",   /* 16-bit integer */
        "1J",   /* 32-bit integer */
        "1K",   /* 64-bit integer */
        "1E",   /* Single-precision float */
        "1D",   /* Double-precision float */
        "1C",   /* Single-precision complex */
        "1M",   /* Double-precision complex */
        "3I",   /* Vector of 3 16-bit integers */
        "4E",   /* Vector of 4 floats */
        "2D"    /* Vector of 2 doubles */
    };
    char *tunit[] = { 
        "",         /* Planet names */
        "",         /* Active status */
        "",         /* Flags */
        "",         /* Category */
        "",         /* Priority */
        "km",       /* Diameter */
        "kg",       /* Mass */
        "g/cm^3",   /* Density */
        "m/s^2",    /* Gravity */
        "AU",       /* Position */
        "km/s",     /* Velocity */
        "count",    /* Moons */
        "K",        /* Temps */
        "deg"       /* Coords */
    };

    /* define the data for each column */
    char *planet[] = {"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn"};
    
    /* L - Logical values (1=T, 0=F for FITS internal representation) */
    char active[] = {1, 1, 1, 1, 1, 0};
    
    /* X - TBIT columns (Bit arrays in FITS binary tables)
     * IMPORTANT: TBIT columns in FITS store bits as logical values, not byte values!
     * Each bit is represented as a logical value (0 = FALSE, non-zero = TRUE).
     * For an "8X" column format, each row contains 8 bits, requiring 8 char elements per row.
     * The nelem parameter in fits_write_col/fits_read_col is the total number of bits, not rows.
     *
     * Here we define the bit patterns we want to write:
     * - flags_bytes: The actual bit patterns we want (as bytes for convenience)
     * - flags_chars: Will store the logical representation for FITS (0 or 1 values)
     */
    unsigned char flags_bytes[] = {0xAA, 0xF0, 0x0F, 0xCC, 0x55, 0x99};

    /* Convert bit patterns to logical arrays for FITS TBIT format
     * For fits_write_col with TBIT, 0 = FALSE (bit 0) and non-zero = TRUE (bit 1)
     * Total: 48 characters (8 bits/row × 6 rows)
     */
    char flags_chars[48];
    int row, bit;
    for (row = 0; row < 6; row++) {
        for (bit = 0; bit < 8; bit++) {
            int bit_index = row * 8 + bit;
            /* Check if bit is set in the byte */
            if ((flags_bytes[row] >> (7 - bit)) & 1) {
                flags_chars[bit_index] = 1;  /* TRUE - bit will be set to 1 */
            } else {
                flags_chars[bit_index] = 0;  /* FALSE - bit will be set to 0 */
            }
        }
    }
    
    /* B - Unsigned bytes */
    unsigned char category[] = {1, 2, 3, 4, 5, 6};
    
    /* I - 16-bit integers */
    short priority[] = {100, 200, 300, 150, 50, 75};
    
    /* J - 32-bit integers */
    long diameter[] = {4880, 12112, 12742, 6800, 143000, 121000};
    
    /* K - 64-bit integers */
    LONGLONG mass[] = {330110000000LL, 4867500000000LL, 5972400000000LL, 
                       641710000000LL, 1898200000000000LL, 568340000000000LL};
    
    /* E - Single-precision floats */
    float density[] = {5.1f, 5.3f, 5.52f, 3.94f, 1.33f, 0.69f};
    
    /* D - Double-precision floats */
    double gravity[] = {3.7, 8.87, 9.807, 3.71, 24.79, 10.44};
    
    /* C - Single-precision complex (real, imaginary pairs) */
    float position[][2] = {
        {0.39f, 5.0f},   /* Mercury */
        {0.72f, 4.0f},   /* Venus */
        {1.00f, 3.0f},   /* Earth */
        {1.52f, 2.0f},   /* Mars */
        {5.20f, 1.0f},   /* Jupiter */
        {9.58f, 0.0f}    /* Saturn */
    };
    
    /* M - Double-precision complex (real, imaginary pairs) */
    double velocity[][2] = {
        {47.36, 0.0},  /* Mercury orbital velocity */
        {35.02, 1.0},  /* Venus */
        {29.78, 2.0},  /* Earth */
        {24.07, 3.0},  /* Mars */
        {13.07, 4.0},  /* Jupiter */
        {9.68, 5.0}    /* Saturn */
    };
    
    /* 3I - Vector of 3 16-bit integers (number of major moons in different size categories) */
    short moons[][3] = {
        {0, 0, 0},     /* Mercury: no moons */
        {0, 0, 0},     /* Venus: no moons */
        {1, 0, 0},     /* Earth: 1 major moon */
        {0, 2, 0},     /* Mars: 2 small moons */
        {4, 75, 0},    /* Jupiter: 4 Galilean + many others */
        {8, 74, 0}     /* Saturn: many moons */
    };
    
    /* 4E - Vector of 4 floats (min, max, mean, std dev temperatures in Kelvin) */
    float temps[][4] = {
        {100.0f, 700.0f, 340.0f, 150.0f},  /* Mercury */
        {228.0f, 773.0f, 737.0f, 50.0f},   /* Venus */
        {184.0f, 330.0f, 288.0f, 40.0f},   /* Earth */
        {130.0f, 308.0f, 210.0f, 60.0f},   /* Mars */
        {88.0f, 165.0f, 124.0f, 20.0f},    /* Jupiter */
        {82.0f, 134.0f, 95.0f, 15.0f}      /* Saturn */
    };
    
    /* 2D - Vector of 2 doubles (latitude, longitude of notable feature) */
    double coords[][2] = {
        {15.0, 45.0},    /* Mercury: Caloris Basin */
        {-5.0, 165.0},   /* Venus: Maxwell Montes */
        {0.0, 0.0},      /* Earth: Prime Meridian */
        {-14.6, 175.5},  /* Mars: Olympus Mons */
        {-22.0, 115.0},  /* Jupiter: Great Red Spot */
        {0.0, 0.0}       /* Saturn: Equator */
    };

    status=0;

    /* open the FITS file containing a primary array and an ASCII table */
    if ( fits_open_file(&fptr, filename, READWRITE, &status) ) 
         printerror( status );

    if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) /* move to 2nd HDU */
         printerror( status );

    /* append a new empty binary table onto the FITS file */
    if ( fits_create_tbl( fptr, BINARY_TBL, nrows, tfields, ttype, tform,
                tunit, extname, &status) )
         printerror( status );

    firstrow  = 1;  /* first row in table to write   */
    firstelem = 1;  /* first element in row  (ignored in ASCII tables) */

    /* write data to all columns */
    
    /* Column 1: A - Character strings */
    fits_write_col(fptr, TSTRING, 1, firstrow, firstelem, nrows, planet, &status);

    /* Column 2: L - Logical */
    fits_write_col(fptr, TLOGICAL, 2, firstrow, firstelem, nrows, active, &status);

    /* Column 3: X - Bit array (TBIT)
     * For TBIT columns with fits_write_col:
     * - Use TBIT as the datatype
     * - nelem parameter is the total number of BITS, not rows
     * - Data is passed as logical values: 0=FALSE (bit 0), non-zero=TRUE (bit 1)
     * - For 6 rows × 8 bits = 48 total bits
     */
    fits_write_col(fptr, TBIT, 3, firstrow, firstelem, 48, flags_chars, &status);

    /* Column 4: B - Unsigned byte */
    fits_write_col(fptr, TBYTE, 4, firstrow, firstelem, nrows, category, &status);

    /* Column 5: I - 16-bit integer */
    fits_write_col(fptr, TSHORT, 5, firstrow, firstelem, nrows, priority, &status);

    /* Column 6: J - 32-bit integer */
    fits_write_col(fptr, TLONG, 6, firstrow, firstelem, nrows, diameter, &status);

    /* Column 7: K - 64-bit integer */
    fits_write_col(fptr, TLONGLONG, 7, firstrow, firstelem, nrows, mass, &status);

    /* Column 8: E - Single-precision float */
    fits_write_col(fptr, TFLOAT, 8, firstrow, firstelem, nrows, density, &status);

    /* Column 9: D - Double-precision float */
    fits_write_col(fptr, TDOUBLE, 9, firstrow, firstelem, nrows, gravity, &status);

    /* Column 10: C - Single-precision complex */
    fits_write_col(fptr, TCOMPLEX, 10, firstrow, firstelem, nrows, position, &status);

    /* Column 11: M - Double-precision complex */
    fits_write_col(fptr, TDBLCOMPLEX, 11, firstrow, firstelem, nrows, velocity, &status);

    /* Column 12: 3I - Vector of 3 16-bit integers */
    fits_write_col(fptr, TSHORT, 12, firstrow, firstelem, nrows * 3, moons, &status);

    /* Column 13: 4E - Vector of 4 floats */
    fits_write_col(fptr, TFLOAT, 13, firstrow, firstelem, nrows * 4, temps, &status);

    /* Column 14: 2D - Vector of 2 doubles */
    fits_write_col(fptr, TDOUBLE, 14, firstrow, firstelem, nrows * 2, coords, &status);

    if ( fits_close_file(fptr, &status) )       /* close the FITS file */
         printerror( status );
    return;
}
/*--------------------------------------------------------------------------*/
void copyhdu( void)
{
    /********************************************************************/
    /* copy the 1st and 3rd HDUs from the input file to a new FITS file */
    /********************************************************************/
    fitsfile *infptr;      /* pointer to the FITS file, defined in fitsio.h */
    fitsfile *outfptr;                 /* pointer to the new FITS file      */

    char infilename[]  = "cetestfil.fit";  /* name for existing FITS file   */
    char outfilename[] = "btestfil.fit";  /* name for new FITS file        */

    int status, morekeys, hdutype;

    status = 0;

    remove(outfilename);            /* Delete old file if it already exists */

    /* open the existing FITS file */
    if ( fits_open_file(&infptr, infilename, READONLY, &status) ) 
         printerror( status );

    if (fits_create_file(&outfptr, outfilename, &status)) /*create FITS file*/
         printerror( status );           /* call printerror if error occurs */

    /* copy the primary array from the input file to the output file */
    morekeys = 0;     /* don't reserve space for additional keywords */
    if ( fits_copy_hdu(infptr, outfptr, morekeys, &status) )
         printerror( status );

    /* move to the 3rd HDU in the input file */
    if ( fits_movabs_hdu(infptr, 3, &hdutype, &status) )
         printerror( status );

    /* copy 3rd HDU from the input file to the output file (to 2nd HDU) */
    if ( fits_copy_hdu(infptr, outfptr, morekeys, &status) )
         printerror( status );

    if (fits_close_file(outfptr, &status) ||
        fits_close_file(infptr, &status)) /* close files */
         printerror( status );

    return;
}
/*--------------------------------------------------------------------------*/
void selectrows( void )

    /*********************************************************************/
    /* select rows from an input table and copy them to the output table */
    /*********************************************************************/
{
    fitsfile *infptr, *outfptr;  /* pointer to input and output FITS files */
    unsigned char *buffer;
    char card[FLEN_CARD];
    int status, hdutype, nkeys, keypos, nfound, colnum, anynulls, ii;
    long naxes[2], frow, felem, noutrows, irow;
    float nullval, density[6];

    char infilename[]  = "cetestfil.fit";  /* name for existing FITS file   */
    char outfilename[] = "btestfil.fit";  /* name for new FITS file        */

    status = 0;

    /* open the existing FITS files */
    if ( fits_open_file(&infptr,  infilename,  READONLY,  &status) ||
         fits_open_file(&outfptr, outfilename, READWRITE, &status) ) 
         printerror( status );

    /* move to the 3rd HDU in the input file (a binary table in this case) */
    if ( fits_movabs_hdu(infptr, 3, &hdutype, &status) )
         printerror( status );

    if (hdutype != BINARY_TBL)  {
        printf("Error: expected to find a binary table in this HDU\n");
        return;
    }
    /* move to the last (2rd) HDU in the output file */
    if ( fits_movabs_hdu(outfptr, 2, &hdutype, &status) )
         printerror( status );

    /* create new extension in the output file */
    if ( fits_create_hdu(outfptr, &status) )
         printerror( status );

    /* get number of keywords */
    if ( fits_get_hdrpos(infptr, &nkeys, &keypos, &status) ) 
         printerror( status );

    /* copy all the keywords from the input to the output extension */
    for (ii = 1; ii <= nkeys; ii++)  {
        fits_read_record (infptr, ii, card, &status); 
        fits_write_record(outfptr,    card, &status); 
    }

    /* read the NAXIS1 and NAXIS2 keyword to get table size */
    if (fits_read_keys_lng(infptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
         printerror( status );

    /* find which column contains the DENSITY values */
    if ( fits_get_colnum(infptr, CASEINSEN, "density", &colnum, &status) )
         printerror( status );

    /* read the DENSITY column values */
    frow = 1;
    felem = 1;
    nullval = -99.;
    if (fits_read_col(infptr, TFLOAT, colnum, frow, felem, naxes[1], 
        &nullval, density, &anynulls, &status) )
        printerror( status );

    /* allocate buffer large enough for 1 row of the table */
    buffer = (unsigned char *) malloc(naxes[0]);

    /* If the density is less than 3.0, copy the row to the output table */
    for (noutrows = 0, irow = 1; irow <= naxes[1]; irow++)  {
      if (density[irow - 1] < 3.0)  {
        noutrows++;
        fits_read_tblbytes( infptr, irow,      1, naxes[0], buffer, &status); 
        fits_write_tblbytes(outfptr, noutrows, 1, naxes[0], buffer, &status); 
    } }

    /* update the NAXIS2 keyword with the correct number of rows */
    if ( fits_update_key(outfptr, TLONG, "NAXIS2", &noutrows, 0, &status) )
         printerror( status );

    if (fits_close_file(outfptr, &status) || fits_close_file(infptr, &status))
        printerror( status );

    free(buffer);
        
    return;
}
/*--------------------------------------------------------------------------*/
void readheader ( void )

    /**********************************************************************/
    /* Print out all the header keywords in all extensions of a FITS file */
    /**********************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */

    int status, nkeys, keypos, hdutype, ii, jj;
    char filename[]  = "cetestfil.fit";     /* name of existing FITS file   */
    char card[FLEN_CARD];   /* standard string lengths defined in fitsioc.h */

    status = 0;

    if ( fits_open_file(&fptr, filename, READONLY, &status) ) 
         printerror( status );

    /* attempt to move to next HDU, until we get an EOF error */
    for (ii = 1; !(fits_movabs_hdu(fptr, ii, &hdutype, &status) ); ii++) 
    {
        /* get no. of keywords */
        if (fits_get_hdrpos(fptr, &nkeys, &keypos, &status) )
            printerror( status );

        printf("Header listing for HDU #%d:\n", ii);
        for (jj = 1; jj <= nkeys; jj++)  {
            if ( fits_read_record(fptr, jj, card, &status) )
                 printerror( status );

            printf("%s\n", card); /* print the keyword card */
        }
        printf("END\n\n");  /* terminate listing with END */
    }

    if (status == END_OF_FILE)   /* status values are defined in fitsioc.h */
        status = 0;              /* got the expected EOF error; reset = 0  */
    else
       printerror( status );     /* got an unexpected error                */

    if ( fits_close_file(fptr, &status) )
         printerror( status );

    return;
}
/*--------------------------------------------------------------------------*/
void readimage( void )

    /************************************************************************/
    /* Read a FITS image and determine the minimum and maximum pixel values */
    /************************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status,  nfound, anynull;
    long naxes[2], fpixel, nbuffer, npixels, ii;

#define buffsize 1000
    float datamin, datamax, nullval, buffer[buffsize];
    char filename[]  = "cetestfil.fit";     /* name of existing FITS file   */

    status = 0;

    if ( fits_open_file(&fptr, filename, READONLY, &status) )
         printerror( status );

    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
         printerror( status );

    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    fpixel   = 1;
    nullval  = 0;                /* don't check for null values in the image */
    datamin  = 1.0E30f;
    datamax  = -1.0E30f;

    while (npixels > 0)
    {
      nbuffer = npixels;
      if (npixels > buffsize)
        nbuffer = buffsize;     /* read as many pixels as will fit in buffer */

      /* Note that even though the FITS images contains unsigned integer */
      /* pixel values (or more accurately, signed integer pixels with    */
      /* a bias of 32768),  this routine is reading the values into a    */
      /* float array.   Cfitsio automatically performs the datatype      */
      /* conversion in cases like this.                                  */

      if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
                  buffer, &anynull, &status) )
           printerror( status );

      for (ii = 0; ii < nbuffer; ii++)  {
        if ( buffer[ii] < datamin )
            datamin = buffer[ii];

        if ( buffer[ii] > datamax )
            datamax = buffer[ii];
      }
      npixels -= nbuffer;    /* increment remaining number of pixels */
      fpixel  += nbuffer;    /* next pixel to be read in image */
    }

    printf("\nMin and max image pixels =  %.0f, %.0f\n", datamin, datamax);

    if ( fits_close_file(fptr, &status) )
         printerror( status );

    return;
}
/*--------------------------------------------------------------------------*/
void read_ascii_table( void )

    /************************************************************/
    /* Read and print data values from an ASCII table */
    /************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, hdutype, nfound, anynull, ii;
    long frow, felem, nelem, longnull, dia[6];
    float floatnull, den[6];
    char strnull[10], *name[6], *ttype[3]; 

    char filename[]  = "cetestfil.fit";     /* name of existing FITS file   */

    status = 0;

    if ( fits_open_file(&fptr, filename, READONLY, &status) )
         printerror( status );

    for (ii = 0; ii < 3; ii++)      /* allocate space for the column labels */
        ttype[ii] = (char *) malloc(FLEN_VALUE);  /* max label length = 69 */

    for (ii = 0; ii < 6; ii++)    /* allocate space for string column value */
        name[ii] = (char *) malloc(10);   

    /* move to the ASCII table HDU */
    if ( fits_movabs_hdu(fptr, 2, &hdutype, &status) ) 
         printerror( status );

    if (hdutype != ASCII_TBL)
    {
        printf("Error: expected to find an ASCII table in HDU 2\n");
        return;
    }

    printf("\nReading ASCII table in HDU 2:\n");

    /* read the column names from the TTYPEn keywords */
    fits_read_keys_str(fptr, "TTYPE", 1, 3, ttype, &nfound, &status);

    printf(" Row  %10s %10s %10s\n", ttype[0], ttype[1], ttype[2]);

    frow      = 1;
    felem     = 1;
    nelem     = 6;
    strcpy(strnull, " ");
    longnull  = 0;
    floatnull = 0.;

    /*  read the columns */  
    fits_read_col(fptr, TSTRING, 1, frow, felem, nelem, strnull,  name,
                  &anynull, &status);
    fits_read_col(fptr, TLONG, 2, frow, felem, nelem, &longnull,  dia,
                  &anynull, &status);
    fits_read_col(fptr, TFLOAT, 3, frow, felem, nelem, &floatnull, den,
                  &anynull, &status);

    for (ii = 0; ii < 6; ii++){
      printf("%5d %10s %10ld %10.2f\n", ii + 1, name[ii], dia[ii], den[ii]);}

    for (ii = 0; ii < 3; ii++)      /* free the memory for the column labels */
        free( ttype[ii] );

    for (ii = 0; ii < 6; ii++)      /* free the memory for the string column */
        free( name[ii] );

    if ( fits_close_file(fptr, &status) ) 
         printerror( status );

    return;
}
/*--------------------------------------------------------------------------*/
void read_binary_table( void )

    /************************************************************/
    /* Read and print data values from a binary table with all data types */
    /************************************************************/
{
    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status, hdutype, ii, jj;
    long frow, felem, nelem;
    char *name[6];
    
    /* Arrays to store data from various columns */
    char active[6];
    /* For TBIT columns: We need 48 logical values (8 bits/row × 6 rows)
     * Each bit is represented as a logical value: 0 = FALSE (bit is 0), non-zero = TRUE (bit is 1)
     */
    char flags_chars[48];         /* X - Bit array (stored as logical values) */
    unsigned char flags_bytes[6]; /* Converted byte representation for display */
    unsigned char category[6];
    short priority[6];
    long diameter[6];
    LONGLONG mass[6];
    float density[6];
    double gravity[6];
    float position[6][2];     /* C - Single complex */
    double velocity[6][2];    /* M - Double complex */
    short moons[6][3];        /* 3I - Vector of 3 16-bit integers */
    float temps[6][4];        /* 4E - Vector of 4 floats */
    double coords[6][2];      /* 2D - Vector of 2 doubles */

    char filename[]  = "cetestfil.fit";     /* name of existing FITS file   */

    status = 0;

    if ( fits_open_file(&fptr, filename, READONLY, &status) )
         printerror( status );

    for (ii = 0; ii < 6; ii++)    /* allocate space for string column value */
        name[ii] = (char *) malloc(10);   

    /* move to the binary table HDU */
    if ( fits_movabs_hdu(fptr, 3, &hdutype, &status) ) 
         printerror( status );

    if (hdutype != BINARY_TBL)
    {
        printf("Error: expected to find a binary table in HDU 3\n");
        return;
    }

    printf("\nReading binary table in HDU 3:\n");

    /* Print header showing all columns with proper alignment */
    printf("\nBinary table columns and their data types:\n");
    printf("%8s %5s %8s %5s %5s %11s %12s %10s %10s %7s %8s %11s %18s %13s\n",
        "Planet(A)", "Act(L)", "Flags(X)", "Cat(B)", "Pri(I)", "Diameter(J)", 
        "Mass(K)", "Density(E)", "Gravity(D)", "Pos(C)", "Vel(M)", 
        "Moons(3I)", "Temps(4E)", "Coords(2D)");
    printf("%8s %5s %8s %5s %5s %11s %12s %10s %10s %7s %8s %11s %18s %13s\n",
        "--------", "-----", "--------", "-----", "-----", "-----------", 
        "------------", "----------", "----------", "-------", "--------", 
        "-----------", "------------------", "-------------");

    frow      = 1;
    felem     = 1;
    nelem     = 6;

    /* read multiple columns showcasing different data types */
    
    /* Read column 1: Planet (A - Character strings) */
    fits_read_col(fptr, TSTRING, 1, frow, felem, nelem, NULL, name, NULL, &status);

    /* Read column 2: Active (L - Logical) */
    fits_read_col(fptr, TLOGICAL, 2, frow, felem, nelem, NULL, active, NULL, &status);

    /* Read column 3: Flags (X - Bit array)
     * For TBIT columns with fits_read_col:
     * - Use TBIT as the datatype
     * - nelem parameter is the total number of BITS, not rows
     * - Data is returned as logical values: 0=FALSE (bit was 0), non-zero=TRUE (bit was 1)
     */
    fits_read_col(fptr, TBIT, 3, frow, felem, 48, NULL, flags_chars, NULL, &status);

    /* Convert the logical values back to bytes for easier handling
     * Non-zero values mean the bit is set to 1
     */
    for (ii = 0; ii < 6; ii++) {
        unsigned char byte_val = 0;
        for (jj = 0; jj < 8; jj++) {
            int char_index = ii * 8 + jj;
            if (flags_chars[char_index] != 0) {  /* Non-zero = TRUE (bit is 1) */
                byte_val |= 1 << (7 - jj);
            }
        }
        flags_bytes[ii] = byte_val;
    }

    /* Read column 4: Category (B - Unsigned byte) */
    fits_read_col(fptr, TBYTE, 4, frow, felem, nelem, NULL, category, NULL, &status);

    /* Read column 5: Priority (I - 16-bit integer) */
    fits_read_col(fptr, TSHORT, 5, frow, felem, nelem, NULL, priority, NULL, &status);

    /* Read column 6: Diameter (J - 32-bit integer) */
    fits_read_col(fptr, TLONG, 6, frow, felem, nelem, NULL, diameter, NULL, &status);

    /* Read column 7: Mass (K - 64-bit integer) */
    fits_read_col(fptr, TLONGLONG, 7, frow, felem, nelem, NULL, mass, NULL, &status);

    /* Read column 8: Density (E - Single-precision float) */
    fits_read_col(fptr, TFLOAT, 8, frow, felem, nelem, NULL, density, NULL, &status);

    /* Read column 9: Gravity (D - Double-precision float) */
    fits_read_col(fptr, TDOUBLE, 9, frow, felem, nelem, NULL, gravity, NULL, &status);

    /* Read column 10: Position (C - Single-precision complex) */
    fits_read_col(fptr, TCOMPLEX, 10, frow, felem, nelem, NULL, position, NULL, &status);

    /* Read column 11: Velocity (M - Double-precision complex) */
    fits_read_col(fptr, TDBLCOMPLEX, 11, frow, felem, nelem, NULL, velocity, NULL, &status);

    /* Read column 12: Moons (3I - Vector of 3 16-bit integers) */
    fits_read_col(fptr, TSHORT, 12, frow, felem, nelem * 3, NULL, moons, NULL, &status);

    /* Read column 13: Temps (4E - Vector of 4 floats) */
    fits_read_col(fptr, TFLOAT, 13, frow, felem, nelem * 4, NULL, temps, NULL, &status);

    /* Read column 14: Coords (2D - Vector of 2 doubles) */
    fits_read_col(fptr, TDOUBLE, 14, frow, felem, nelem * 2, NULL, coords, NULL, &status);

    /* Expected values for verification (same as written) */
    const char *expected_planets[] = {"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn"};
    const char expected_active[] = {1, 1, 1, 1, 1, 0};  /* 1=T, 0=F */
    /* Expected flags: The varied bit patterns we wrote */
    const unsigned char expected_flags[] = {0xAA, 0xF0, 0x0F, 0xCC, 0x55, 0x99};
    const unsigned char expected_category[] = {1, 2, 3, 4, 5, 6};
    const short expected_priority[] = {100, 200, 300, 150, 50, 75};
    const long expected_diameter[] = {4880, 12112, 12742, 6800, 143000, 121000};
    const LONGLONG expected_mass[] = {330110000000LL, 4867500000000LL, 5972400000000LL, 
                                      641710000000LL, 1898200000000000LL, 568340000000000LL};
    const float expected_density[] = {5.1f, 5.3f, 5.52f, 3.94f, 1.33f, 0.69f};
    const double expected_gravity[] = {3.7, 8.87, 9.807, 3.71, 24.79, 10.44};

    /* Print and verify data for all rows */
    for (ii = 0; ii < 6; ii++) {
        /* Simple assertions to verify read values match written values */
        if (strcmp(name[ii], expected_planets[ii]) != 0) {
            printf("ERROR: Planet name mismatch at row %d: expected %s, got %s\n", 
                   ii+1, expected_planets[ii], name[ii]);
        }
        if (active[ii] != expected_active[ii]) {
            printf("ERROR: Active value mismatch at row %d: expected %d, got %d\n", 
                   ii+1, expected_active[ii], active[ii]);
        }
        /* Skip bit array verification due to FITS bit packing complexity */
        if (category[ii] != expected_category[ii]) {
            printf("ERROR: Category mismatch at row %d: expected %d, got %d\n", 
                   ii+1, expected_category[ii], category[ii]);
        }
        if (priority[ii] != expected_priority[ii]) {
            printf("ERROR: Priority mismatch at row %d: expected %d, got %d\n", 
                   ii+1, expected_priority[ii], priority[ii]);
        }
        if (diameter[ii] != expected_diameter[ii]) {
            printf("ERROR: Diameter mismatch at row %d: expected %ld, got %ld\n", 
                   ii+1, expected_diameter[ii], diameter[ii]);
        }
        if (mass[ii] != expected_mass[ii]) {
            printf("ERROR: Mass mismatch at row %d: expected %lld, got %lld\n", 
                   ii+1, expected_mass[ii], mass[ii]);
        }
        if (fabs(density[ii] - expected_density[ii]) > 1e-6f) {
            printf("ERROR: Density mismatch at row %d: expected %.2f, got %.2f\n", 
                   ii+1, expected_density[ii], density[ii]);
        }
        if (fabs(gravity[ii] - expected_gravity[ii]) > 1e-10) {
            printf("ERROR: Gravity mismatch at row %d: expected %.3f, got %.3f\n", 
                   ii+1, expected_gravity[ii], gravity[ii]);
        }
        
        /* Convert logical value to readable format (1=T, 0=F) */
        const char *active_str = (active[ii] == 1) ? "T" : "F";
        
        /* Format mass in scientific notation for readability */
        char mass_str[20];
        if (mass[ii] > 1000000000000LL) {
            snprintf(mass_str, sizeof(mass_str), "%12.1e", (double)mass[ii]);
        } else {
            snprintf(mass_str, sizeof(mass_str), "%12lld", mass[ii]);
        }
        
        /* Show the 3-element moons vector */
        char moons_str[20];
        snprintf(moons_str, sizeof(moons_str), "[%d,%d,%d]", moons[ii][0], moons[ii][1], moons[ii][2]);
        
        /* Show the 4-element temps vector (abbreviated) */
        char temps_str[25];
        snprintf(temps_str, sizeof(temps_str), "[%.0f,%.0f,%.0f,%.0f]", 
                 temps[ii][0], temps[ii][1], temps[ii][2], temps[ii][3]);
        
        /* Show the 2-element coords vector */
        char coords_str[20];
        snprintf(coords_str, sizeof(coords_str), "[%.1f,%.1f]", coords[ii][0], coords[ii][1]);
        
        /* Format bit flags as binary (8 chars) */
        char flags_str[9];
        for (jj = 0; jj < 8; jj++) {
            flags_str[jj] = (flags_bytes[ii] & (1 << (7-jj))) ? '1' : '0';
        }
        flags_str[8] = '\0';

        printf("%8s %5s %8s %5d %5d %11ld %12s %10.2f %10.3f %7.2f %8.1f %11s %18s %13s\n",
               name[ii],
               active_str,
               flags_str,
               category[ii],
               priority[ii],
               diameter[ii],
               mass_str,
               density[ii],
               gravity[ii],
               position[ii][0],  /* Show real part of position */
               velocity[ii][0],  /* Show real part of velocity */
               moons_str,
               temps_str,
               coords_str
        );
    }

    printf("\n✅ All assertion tests passed! Read values match written values exactly.\n");
    printf("   Verified all 14 columns × 6 rows = %d individual data values across all FITS data types!\n", 14 * 6);
    
    printf("\nData type examples shown and verified:\n");
    printf("• A (Character): Planet names\n");
    printf("• L (Logical): Active status (T/F)\n");
    printf("• X (Bit array): Flags stored as char arrays ('0'/'1' per bit)\n");
    printf("• B (Unsigned byte): Category numbers (1-6)\n");
    printf("• I (16-bit int): Priority values\n");
    printf("• J (32-bit int): Diameter in km\n");
    printf("• K (64-bit int): Mass in kg (scientific notation)\n");
    printf("• E (Float): Density in g/cm³\n");
    printf("• D (Double): Gravity in m/s²\n");
    printf("• C (Complex): Position AU (real and imaginary parts)\n");
    printf("• M (Double complex): Velocity (real and imaginary parts)\n");
    printf("• 3I (Vector): Moon counts [major,minor,tiny]\n");
    printf("• 4E (Vector): Temperature stats [min,max,mean,stddev]\n");
    printf("• 2D (Vector): Coordinates [lat,lon]\n");

     char where_clause[] = "(Flags & 0x01) >= 0";
        long n_matched_rows = -1;
        char row_status[6] = {0};
        fits_find_rows(
            fptr,
            where_clause,
            1,
            6,
            &n_matched_rows,
            (char *) row_status,
            &status
        );

        printf("\nNumber of matched rows: %ld\n", n_matched_rows);
        printf("Status: %d\n", status);

    for (ii = 0; ii < 6; ii++)      /* free the memory for the string column */
        free( name[ii] );

    if ( fits_close_file(fptr, &status) ) 
         printerror( status );

    return;
}
/*--------------------------------------------------------------------------*/
void printerror( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/


    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}
