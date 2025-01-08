# How this was converted into a Calkit project

Inside the project directory:

```sh
calkit new project . --title "Simulating a NACA airfoil with OpenFOAM" --no-commit
```

```sh
calkit new docker-env \
  --name foam \
  --from microfluidica/openfoam:2406 \
  --add-layer miniforge \
  --add-layer foampy \
  --image openfoam-v2406-foampy
```
