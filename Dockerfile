FROM python:3.6

ADD . /app
WORKDIR /app

RUN apt-get -qq update && apt-get install -y \
         gfortran \
         gnuplot \
         libblas-dev \
         liblapack-dev \
         libplot-dev \
         plotutils \
    && rm -rf /var/lib/apt/lists/*

RUN pip install -r requirements.txt

## biggles compilation cannot see numpy when invoked in same pip command
RUN pip install biggles

ENV PYTHONPATH $PYTHONPATH:`pwd`

## install TM-Align
RUN cd /tmp && wget http://zhanglab.ccmb.med.umich.edu/TM-align/TMalign.gz \
    && gzip -d TMalign.gz \
    && chmod +x TMalign \
    && mv TMalign /usr/local/bin/

## install DSSP, Pymol, surfrace dependencies
RUN apt-get -qq update && apt-get install -y \
         dssp \
         pymol \
         unzip \  # required for surfrace installation
         libstdc++5 \ # old libc version required for surfrace binary
    && rm -rf /var/lib/apt/lists/*

## install SurfaceRacer
RUN cd /tmp \
    && wget http://pharmacy.uky.edu/sites/pharmacy.uky.edu/files/files/tsodikov/surface_racer_5.0_linux64.zip \
    && unzip surface_racer_5.0_linux64.zip \
    && rm surface_racer_5.0_linux64.zip \
    && mv surface_racer_5.0_64bit/surfrace5_0_linux_64bit /usr/local/bin/surfrace \
    && chmod +x /usr/local/bin/surfrace \
    && rm -rf surface_racer*

CMD ["python", "biskit/test.py", "-e", "exe"]

