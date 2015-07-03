Compiler

gcc-5 -fopenmp -O3 viho_alt.c -lfftw3 -lX11 -ltiff -ljpeg -lpng

compile with openMP

dans la version alternative, non interactive, on donne l'image puis les 9 coeffs ; viho_alt sauvegarde le résultat dans une nouvelle image (au format .png, de taille WOUT*HOUT = 512*512)

dans la version démo, non interactive, on donne l'image, puis quatre points (au format abscisse ordonnée) puis les quatre points sur lesquels ont souhaite les envoyer ; viho_demo sauvegarde le résultat dans une nouvelle image (au format .png, de taille WOUT*HOUT = 512*512)
exemple : tissu_a est un tapis dont les sommets, listés dans le sens horaire, ont pour coordonnées (495,5), (898,77), (596,634) et (10,261) donc l'appel
  ./a tissu_a.jpg 495 5 898 77 569 634 10 261 0 0 512 0 512 512 0 512
  redresse le tapis en un carré 512*512
