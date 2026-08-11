/*
 * Differential oracle: runs tests/fixtures/eval_corpus.txt through the real
 * CFITSIO expression parser and prints the same one-line summary that
 * tests/test_eval_corpus.rs prints, so the output can be diffed directly
 * against tests/fixtures/eval_corpus.golden.
 *
 * See README.md in this directory for how to build and run it.
 *
 * Everything below must stay in step with `create_corpus_table` and `probe`
 * in tests/test_eval_corpus.rs: same columns, same values, same keywords,
 * same rendering.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "fitsio.h"

#define NROWS    3   /* rows in the corpus table                        */
#define MAX_ELEM 8   /* cap on vector elements evaluated per expression */

static int is_integral(int dt)
{
	return dt == TLOGICAL || dt == TBYTE || dt == TSHORT
	    || dt == TINT     || dt == TLONG || dt == TLONGLONG;
}

/* Expressions whose value is not reproducible run parse-only. */
static int is_nondeterministic(const char *expr)
{
	size_t i;
	char up[4096];
	for (i = 0; expr[i] && i < sizeof up - 1; i += 1)
		up[i] = (char)toupper((unsigned char)expr[i]);
	up[i] = 0;
	return strstr(up, "RANDOM") != NULL;
}

static fitsfile *create_corpus_table(void)
{
	fitsfile *f;
	int status = 0;
	char *ttype[] = { "INTCOL", "FLOATCOL", "STRCOL", "BOOLCOL",
			  "DBLCOL", "VECCOL", "BITS", "MATRIX" };
	char *tform[] = { "1J", "1E", "10A", "1L", "1D", "3E", "8X", "6E" };
	long   intcol[]   = { 7, -3, 10 };
	float  floatcol[] = { 2.5f, 4.0f, 0.5f };
	char  *strcol[]   = { "abc", "de", "fghij" };
	char   boolcol[]  = { 1, 0, 1 };
	double dblcol[]   = { 1.25, -2.5, 8.0 };
	float  veccol[]   = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	unsigned char bits[] = { 0xF0, 0x0F, 0xAA };
	float  matrix[]   = { 1, 2, 3, 4, 5, 6, 11, 12, 13,
			      14, 15, 16, 21, 22, 23, 24, 25, 26 };
	long   ikey = 42;
	double dkey = 2.5;
	int    lkey = 1;
	int i;

	fits_create_file(&f, "mem://corpus.fits", &status);
	fits_write_imghdr(f, BYTE_IMG, 0, NULL, &status);
	fits_create_tbl(f, BINARY_TBL, NROWS, 8, ttype, tform, NULL, NULL, &status);

	fits_write_key_str(f, "TDIM8", "(2,3)", NULL, &status);
	fits_write_key_lng(f, "INTKEY", ikey, NULL, &status);
	fits_write_key_dbl(f, "DBLKEY", dkey, 4, NULL, &status);
	fits_write_key_log(f, "LOGKEY", lkey, NULL, &status);
	fits_write_key_str(f, "STRKEY", "hello", NULL, &status);

	fits_write_col(f, TLONG,    1, 1, 1, 3, intcol,   &status);
	fits_write_col(f, TFLOAT,   2, 1, 1, 3, floatcol, &status);
	for (i = 0; i < 3; i += 1)
		fits_write_col(f, TSTRING, 3, i + 1, 1, 1, &strcol[i], &status);
	fits_write_col(f, TLOGICAL, 4, 1, 1, 3, boolcol,  &status);
	fits_write_col(f, TDOUBLE,  5, 1, 1, 3, dblcol,   &status);
	fits_write_col(f, TFLOAT,   6, 1, 1, 9, veccol,   &status);
	fits_write_col(f, TBYTE,    7, 1, 1, 3, bits,     &status);
	fits_write_col(f, TFLOAT,   8, 1, 1, 18, matrix,  &status);

	if (status) {
		fprintf(stderr, "corpus table setup failed: %d\n", status);
		exit(1);
	}
	return f;
}

static void probe(fitsfile *f, const char *expr, char *out, size_t outsz)
{
	int dt = 0, naxis = 0, status = 0, anynul = 0;
	long nelem = 0, naxes[5] = { 0 };
	long n, i;
	char *p = out;
	char *end = out + outsz;

	fits_test_expr(f, (char *)expr, 5, &dt, &nelem, &naxis, naxes, &status);
	if (status) {
		char m[81];
		p += snprintf(p, (size_t)(end - p), "ERR %d", status);
		/* ORACLE_MSGS=1 appends CFITSIO's own error text, which is what
		   a library consumer reads back with fits_read_errmsg. Off by
		   default so the output still matches the golden file. */
		while (fits_read_errmsg(m))
			if (getenv("ORACLE_MSGS"))
				p += snprintf(p, (size_t)(end - p), " |%s", m);
		fits_clear_errmsg();
		return;
	}

	p += snprintf(p, (size_t)(end - p), "OK dt=%d nelem=%ld naxis=%d",
		      dt, nelem, naxis);
	if (naxis > 0 && naxis <= 5) {
		p += snprintf(p, (size_t)(end - p), " naxes=[");
		for (i = 0; i < naxis; i += 1)
			p += snprintf(p, (size_t)(end - p), "%s%ld",
				      i ? ", " : "", naxes[i]);
		p += snprintf(p, (size_t)(end - p), "]");
	}

	/* String results carry no numeric value; the corpus exercises string
	   semantics through the boolean comparison lines instead. */
	if (dt == TSTRING || is_nondeterministic(expr))
		return;

	n = nelem < 1 ? 1 : (nelem > MAX_ELEM ? MAX_ELEM : nelem);

	if (is_integral(dt)) {
		LONGLONG *v = calloc((size_t)(n * NROWS), sizeof *v);
		fits_calc_rows(f, TLONGLONG, (char *)expr, 1, n * NROWS,
			       NULL, v, &anynul, &status);
		if (status) {
			snprintf(p, (size_t)(end - p), " | EVALERR %d", status);
			fits_clear_errmsg();
			free(v);
			return;
		}
		p += snprintf(p, (size_t)(end - p), " | [");
		for (i = 0; i < n * NROWS; i += 1)
			p += snprintf(p, (size_t)(end - p), "%s%lld",
				      i ? ", " : "", v[i]);
		p += snprintf(p, (size_t)(end - p), "]");
		free(v);
	} else {
		double *v = calloc((size_t)(n * NROWS), sizeof *v);
		fits_calc_rows(f, TDOUBLE, (char *)expr, 1, n * NROWS,
			       NULL, v, &anynul, &status);
		if (status) {
			snprintf(p, (size_t)(end - p), " | EVALERR %d", status);
			fits_clear_errmsg();
			free(v);
			return;
		}
		p += snprintf(p, (size_t)(end - p), " | [");
		for (i = 0; i < n * NROWS; i += 1)
			p += snprintf(p, (size_t)(end - p), "%s%.6f",
				      i ? ", " : "", v[i]);
		p += snprintf(p, (size_t)(end - p), "]");
		free(v);
	}

	if (anynul)
		snprintf(p, (size_t)(end - p), " (null)");
}

/* A corpus line is a comment when it is "#" alone or starts with "# ".
   Bare "#NAME" lines are expressions (the #KEYWORD constant form). */
static int is_comment(const char *l)
{
	return strcmp(l, "#") == 0 || strncmp(l, "# ", 2) == 0;
}

int main(int argc, char **argv)
{
	fitsfile *f;
	FILE *in = stdin;
	char line[4096], out[8192];
	int status = 0;

	if (argc > 1 && !(in = fopen(argv[1], "r"))) {
		perror(argv[1]);
		return 1;
	}

	f = create_corpus_table();
	while (fgets(line, sizeof line, in)) {
		size_t n = strcspn(line, "\n");
		line[n] = 0;
		while (n > 0 && (line[n - 1] == ' ' || line[n - 1] == '\t'))
			line[--n] = 0;
		if (!line[0] || is_comment(line))
			continue;
		probe(f, line, out, sizeof out);
		printf("%s\t%s\n", line, out);
	}
	fits_close_file(f, &status);
	return 0;
}
