#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static char *disc[]={
  "DISCLAIMER: The authors do not guarantee that the "
  "implementation of this recipe leads to an eatable dish, for it is "
  "clearly \"system dependent\".",
  NULL
};

static char *es[]={
  "PULPO A FEIRA:\n\nIngredientes: Para 4 personas\n"
	"  - 2 kg. de pulpo\n"
	"  - 2 kg. de patatas\n"
  "  - 100 grs. de pimentón picante\n"
  "  - 100 grs. de sal gorda\n"
  "  - aceite\n\n"
	"Preparación: Se lava el pulpo en agua fría, se pone una olla de cobre "
	"con agua al fuego y cuando rompa a hervir se coge el pulpo, se mete y "
	"se saca del agua tres veces dejando que en cada intervalo vuelva a "
	"hervir el agua. Se deja cocer el pulpo durante unos 20 minutos "
	"retirándolo del fuego y dejándolo reposar durante 5 minutos. A "
	"continuación, se quita del agua y se corta en trozos finos con "
	"unas tijeras. Para servirlo se pone en unos platos de madera "
	"condimentándolo por este orden: sal, pimentón, aceite y se añaden "
	"unos cachelos (Patatas).",
  NULL
};

static char *en[]={
	"OCTOPUS a la GALLEGA:\n\nIngredients: For 4 persons\n"
	"  - 2 kg. of octopus\n"
	"  - 2 kg. of potatoes\n"
  "  - 100 grs. of paprika\n"
  "  - 100 grs. thick salt\n"
  "  - olive oil\n\n"
	"Preparation: Wash the octopus in cold water. Place a cooper pan "
	"with water to the fire, and when it starts boiling submerge the "
	"octopus, and remove it again. Repeat this procedure 3 times, "
	"always waiting for the water to start boiling. Then, boil the "
	"octopus for 20 minutes, remove the pan from the fire and let it "
	"rest for 5 minutes. Remove the octopus from the water, and cut it "
	"in thin slices with some scissors. It is served in wood plates, "
	"spiced in the following order: first the salt, then the paprika "
	"and the olive oil. It is served with potatoes.",
	NULL
};

#define NCOLS 75
#define NO_LANGS 2
static char **rec[] = {disc,en,es};

void F90_FUNC(oct_printrecipe, OCT_PRINTRECIPE)
		 (int *lang)
{
	int i, j, n;
	char *s, c[NCOLS+5];

	/* initialize random numbers */
	srand((unsigned int)time(NULL));

	if(*lang<0 || *lang>=NO_LANGS){
		fprintf(stderr, "Unsupported language (%d)", *lang);
		return;
	}

	/* count recipes in this language */
	n = 0;
	while(rec[*lang][n] != NULL) n++;
	i = (int) (((float)n)*rand()/(RAND_MAX+1.0));

	/* we now print it in 80 column format */
	s = rec[*lang][i];
	do{
		/* skip inital white space */
		for(; *s!='\0' && *s==' '; s++);
		for(i=0; i<NCOLS+1 && *s!='\0' && *s!='\n'; i++)
			c[i] = *s++;
		if(i==NCOLS+1){
			for(; i!=0 && c[i] != ' '; i--, s--);
			if(c[i] == ' ') i--;
		}else{
			c[i] = *s++;
		}
		c[i+1] = '\0';

		if(c[i]!='\0' && c[i] != '\n'){
			n = NCOLS-strlen(c);
			/* add extra spaces */
			for(j=NCOLS-1; j>=0 && i>=0 && n>0; j--, i--){
				c[j] = c[i];
				if(c[j] == ' '){
					c[--j] = ' ';
					n--;
				}
			}
			if(n > 0) /* unlikely, we do not care so much */
				for(i=0; i<n; i++)
					c[i] = ' ';
			c[NCOLS]='\n'; c[NCOLS+1] = '\0';
		}
		printf("%s", c);
	}while(*(s-1));
	printf("\n\n");
	fflush(stdout);
}
