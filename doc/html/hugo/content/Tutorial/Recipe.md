---
title: "Recipe"
tags: ["Basic", "Recipe"]
series: "Tutorial"
---


If the input file sets

```text
 {{< Variable2 "CalculationMode" >}} = recipe
```

then the code will print you a tasty recipe, randomly selected from those available in the package. For example:

```text
 ************************** Calculation Mode **************************
 Input: [CalculationMode = recipe]
 **********************************************************************
 
 Info: Octopus initialization completed.
 Info: Starting calculation mode.
 
 
 OCTOPUS a la GALLEGA:
 
 Ingredients: For 4 persons
   - 2 kg. of octopus
   - 2 kg. of potatoes
   - 100 grs. of paprika
   - 100 grs. thick salt
   - olive oil
 
 Preparation: Wash the octopus in cold water. Place a copper (Cu) pan with water
 on the fire, and when it starts boiling submerge the octopus, and remove
 it again. Repeat this procedure 3 times, always waiting for the water to
 start boiling. Then, boil the octopus for 20 minutes, remove the pan from
 the fire and let it rest for 5 minutes. Remove the octopus from the water,
 and cut it in thin slices with some scissors. It is served on wood plates,
 spiced in the following order: first the salt, then the paprika and the
 olive oil. It is served with potatoes.
 
 
 DISCLAIMER: The authors do not guarantee that the implementation of this
 recipe leads to an edible dish, for it is clearly "system-dependent".
```

There are currently recipes available in English, Spanish, Italian, and Basque. Which one you get is controlled by the LANG environment variable, so if you run "LANG=es octopus", you may get:

```text
 PULPO A FEIRA:
 
 Ingredientes: Para 4 personas
   - 2 kg. de pulpo
   - 2 kg. de papas / patatas
   - 100 grs. de ají / chile / pimentón picante
   - 100 grs. de sal gorda
   - aceite
 
 Preparación: Se lava el pulpo en agua fría, se pone una olla de cobre con
 agua al fuego y cuando rompa a hervir se coge el pulpo, se mete y se saca
 del agua tres veces dejando que en cada intervalo vuelva a hervir el agua.
 Se deja cocer el pulpo durante unos 20 minutos retirándolo del fuego y
 dejándolo reposar durante 5 minutos. A continuación, se quita del agua y
 se corta en trozos finos con unas tijeras. Para servirlo se pone en unos
 platos de madera condimentándolo por este orden: sal, pimentón, aceite y
 se añaden unos cachelos (Patatas).
```

Contributions in other languages welcome!

{{Tutorial_foot|series=Octopus basics|prev=Time-dependent propagation|next=}}



---------------------------------------------
