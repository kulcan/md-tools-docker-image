# MDTools para Docker
Especificación de una imagen en Docker sobre [Ubuntu 22.04 LTS](https://ubuntu.com/download) con herramientas para dinámica molecular. Tamaño aproximado 9.34GB con todos los paquetes.

Esta imagen incluye:
* [Psi4](https://psicode.org/)
* [ASE](https://wiki.fysik.dtu.dk/ase/index.html)
* [Sympy](https://www.sympy.org/)
* [OpenMM](https://openmm.org/)

Todo corriendo bajo un mismo entorno en [Python 3.8](https://www.python.org/) y Mini Conda [Conda](https://docs.conda.io/en/latest/).

## Instrucciones de instalación
1. Tener [Docker](https://www.docker.com/) instalado en su computadora.

2. Para crear la imagen simplemente hay que ejecutar el comando

    ```
    docker build . -t md-tools:v11 --progress=plain
    ```

3. Despúes podemos correr el contenedor en el puerto 8888 de nuestra computadora

    ```
    docker run -d -p 8888:8888 md-tools:v11
    ```
    
4. Hecho esto dirigase a [localhost:8888](http://localhost:8888) en su navegador y usted será capaz de ver el jupyterlab del contenedor.

5. Para comprobar el correcto funcionamiento del contenedor ejecute el notebook **test.ipynb**.
