FROM python:3.6

RUN apt-get -qq update && apt-get install -y \
         gfortran \
         gnuplot \
         libblas-dev \
         liblapack-dev \
         libplot-dev \
         plotutils \
    && rm -rf /var/lib/apt/lists/*

## don't use -r requirements.txt here because we have not ADDed file yet
RUN pip install numpy scipy biopython

## biggles compilation cannot see numpy when invoked in same pip command
RUN pip install biggles

## install TM-Align
RUN cd /tmp && wget http://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.gz \
    && gzip -d TMalign.gz \
    && chmod +x TMalign \
    && mv TMalign /usr/local/bin/

## install DSSP, Pymol, surfrace dependencies
RUN apt-get -qq update && apt-get install -y \
         dssp \
         pymol \
         unzip \ 
         libstdc++5 \ 
    && rm -rf /var/lib/apt/lists/*

## install SurfaceRacer
RUN cd /tmp \
    && wget http://pharmacy.uky.edu/sites/pharmacy.uky.edu/files/files/tsodikov/surface_racer_5.0_linux64.zip \
    && unzip surface_racer_5.0_linux64.zip \
    && rm surface_racer_5.0_linux64.zip \
    && mv surface_racer_5.0_64bit/surfrace5_0_linux_64bit /usr/local/bin/surfrace \
    && chmod +x /usr/local/bin/surfrace \
    && rm -rf surface_racer*

## install reduce
RUN cd /tmp \
    && wget "http://kinemage.biochem.duke.edu/php/downlode.php?filename=/downloads/software/reduce31/reduce.3.23.130521.linuxi386.gz" \
    && mv downlode*gz reduce.gz \
    && gzip -d reduce.gz \
    && chmod +x reduce \
    && mv reduce /usr/local/bin/ 

## Everything following ADD is not cached by Docker
ADD . /app
WORKDIR /app

## duplicate, just in case requirements was updated without updating Dockerfile
RUN pip install -r requirements_extended.txt

ENV PYTHONPATH $PYTHONPATH:`pwd`

CMD ["python", "biskit/test.py"]

