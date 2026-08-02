/*
 * Token-level oracle: the pattern definitions of CFITSIO's eval.l with every
 * action replaced by a print, so the token stream can be compared against
 * src/parser/lexer.rs without needing a FITS file or a ParseData.
 *
 * The definitions block below is copied verbatim from ~/code/cfitsio/eval.l.
 * Regenerate this file if eval.l changes.  Rule ORDER is significant: flex
 * takes the longest match and breaks ties by declaration order.
 *
 *   make lexdump && echo 'b101 1.5e3 T TT' | ./lexdump
 */
%{
#include <stdio.h>
%}
%option noyywrap
bit		([bB][01xX]+)
oct		([oO][01234567xX]+)
hex		([hH][0123456789aAbBcCdDeEfFxX]+)
bitconst        (0b[01]+)
hexconst        (0x[0123456789aAbBcCdDeEfF]+)
octconst        (0o[01234567]+)
integer		[0-9]+
boolean         (t|f|T|F)
real		([0-9]*"."[0-9]+)|([0-9]*"."*[0-9]+[eEdD][+-]?[0-9]+)|([0-9]+".")
constant        ("#"[a-zA-Z0-9_]+)|("#""$"[^\n]*"$")
string		([\"][^\"\n]*[\"])|([\'][^\'\n]*[\'])
variable	([a-zA-Z_][a-zA-Z0-9_]*)|("$"[^$\n]*"$")
function	[a-zA-Z][a-zA-Z0-9]+"("
intcast		("(int)"|"(INT)")
fltcast		("(float)"|"(FLOAT)"|"(double)"|"(DOUBLE)")
power		("^"|"**")
not             ("!"|".not."|".NOT."|"not."|"NOT.")
or              ("||"|".or."|".OR."|"or."|"OR.")
and             ("&&"|".and."|".AND."|"and."|"AND.")
equal		("=="|".eq."|".EQ."|"eq."|"EQ.")
not_equal	("!="|".ne."|".NE."|"ne."|"NE.")
greater         (">"|".gt."|".GT."|"gt."|"GT.")
lesser          ("<"|".lt."|".LT."|"lt."|"LT.")
greater_eq	(">="|"=>"|".ge."|".GE."|"ge."|"GE.")
lesser_eq	("<="|"=<"|".le."|".LE."|"le."|"LE.")
xor             ("^^"|".xor."|".XOR.")
nl		\n
%%
[ \t]+     ;
{bit}           { printf("BIT<%s> ", yytext); }
{oct}           { printf("OCT<%s> ", yytext); }
{hex}           { printf("HEX<%s> ", yytext); }
{bitconst}      { printf("BITCONST<%s> ", yytext); }
{octconst}      { printf("OCTCONST<%s> ", yytext); }
{hexconst}      { printf("HEXCONST<%s> ", yytext); }
{integer}       { printf("LONG<%s> ", yytext); }
{boolean}       { printf("BOOLEAN<%s> ", yytext); }
{real}          { printf("DOUBLE<%s> ", yytext); }
{constant}      { printf("CONST<%s> ", yytext); }
{string}        { printf("STRING<%s> ", yytext); }
{variable}      { printf("VAR<%s> ", yytext); }
{function}      { printf("FUNC<%s> ", yytext); }
{intcast}       { printf("INTCAST "); }
{fltcast}       { printf("FLTCAST "); }
{power}         { printf("POWER "); }
{not}           { printf("NOT "); }
{or}            { printf("OR "); }
{and}           { printf("AND "); }
{equal}         { printf("EQ "); }
{not_equal}     { printf("NE "); }
{greater}       { printf("GT "); }
{lesser}        { printf("LT "); }
{greater_eq}    { printf("GTE "); }
{lesser_eq}     { printf("LTE "); }
{xor}           { printf("XOR "); }
{nl}          { printf("\n"); }
.             { printf("CH<%s> ", yytext); }
%%
int main(void) { yylex(); return 0; }
