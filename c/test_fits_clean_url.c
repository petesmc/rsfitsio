#include <stdio.h>
#include <string.h>
#include "fitsio2.h"

// Mock implementation of fits_clean_url for testing
int fits_clean_url(char *input, char *output, int *status);

void test_fits_clean_url() {
    char output[FLEN_FILENAME];
    int status = 0;

    printf("Testing fits_clean_url function...\n");
    // Initialize the output buffer
    memset(output, 0, sizeof(output));


    // Test case 1: Normalizing a simple relative path
    printf("Test case 1: Normalizing a simple relative path\n");
    char input1[] = "/dir/../file.fits";
    fits_clean_url(input1, output, &status);
    if (strcmp(output, "/file.fits") == 0 && status == 0) {
        printf("Test 1 passed.\n");
    } else {
        printf("Test 1 failed. Output: %s, Status: %d\n", output, status);
    }

    // Test case 2: Normalizing an absolute path
    printf("Test case 2: Normalizing an absolute path\n");
    char input2[] = "/home/user/../file.fits";
    fits_clean_url(input2, output, &status);
    if (strcmp(output, "/home/file.fits") == 0 && status == 0) {
        printf("Test 2 passed.\n");
    } else {
        printf("Test 2 failed. Output: %s, Status: %d\n", output, status);
    }

    // Test case 3: Handling redundant slashes
    printf("Test case 3: Handling redundant slashes\n");
    char input3[] = "/home//user///file.fits";
    fits_clean_url(input3, output, &status);
    if (strcmp(output, "/home/user/file.fits") == 0 && status == 0) {
        printf("Test 3 passed.\n");
    } else {
        printf("Test 3 failed. Output: %s, Status: %d\n", output, status);
    }

    // Test case 4: Empty input
    printf("Test case 4: Empty input\n");
    char input4[] = "";
    fits_clean_url(input4, output, &status);
    if (strcmp(output, "") == 0 && status == 0) {
        printf("Test 4 passed.\n");
    } else {
        printf("Test 4 failed. Output: %s, Status: %d\n", output, status);
    }

    // Test case 5: Input with only dots
    printf("Test case 5: Input with only dots\n");
    char input5[] = "././.";
    fits_clean_url(input5, output, &status);
    if (strcmp(output, "") == 0 && status == 0) {
        printf("Test 5 passed.\n");
    } else {
        printf("Test 5 failed. Output: %s, Status: %d\n", output, status);
    }

    // Test case 6: Input with invalid characters
    printf("Test case 6: Input with invalid characters\n");
    char input6[] = "dir/../file?.fits";
    fits_clean_url(input6, output, &status);
    if (strcmp(output, "file?.fits") == 0 && status == 0) {
        printf("Test 6 passed.\n");
    } else {
        printf("Test 6 failed. Output: %s, Status: %d\n", output, status);
    }

    }

int main() {
    test_fits_clean_url();
    return 0;
}
