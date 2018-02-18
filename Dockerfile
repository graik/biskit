FROM python:3.6

## dependencies for Biskit and Biggles
RUN apt-get -qq update && apt-get install -y \
         gfortran \
         gnuplot \
         libblas-dev \
         liblapack-dev \
         libplot-dev \
         plotutils \
    && rm -rf /var/lib/apt/lists/*

## don't use -r requirements.txt here because we have not ADDed file yet
RUN pip install numpy scipy biopython \
## biggles compilation cannot see numpy when invoked in same pip command
    && pip install biggles

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

## install AmberTools 17
RUN cd /tmp \
    && wget http://ambermd.org/downloads/install_ambertools.sh \
    && bash install_ambertools.sh --prefix /opt --non-conda \
    && rm -r /opt/ambertools*bz2 

ENV AMBERHOME /opt/amber17
ENV PATH $PATH:$AMBERHOME/bin

## Everything following ADD is not cached by Docker
## -> the following RUN commands have to be re-executed with every build
ADD . /app
WORKDIR /app

## Now install programs that require registration or cannot be automatically
## downloaded for other reasons
## download manually and put copies into the `downloads` folder

## Try installing Delphi from downloads copy
RUN if test -e downloads/DelPhi_Linux_SP_F95.tar.gz; then \
       mv downloads/DelPhi_Linux_SP_F95.tar.gz /tmp; \
       cd /tmp ; \
       tar xvfz DelPhi_Linux_SP_F95.tar.gz ; \
       mv DelPhi*/ /opt/delphi_sp ; \
       ln -s /opt/delphi_sp/executable/delphi95 /usr/local/bin/delphi ; \
       cd /app ; \
       rm /tmp/DelPhi_Linux_* ; \
       echo "Delphi installed from downloads copy."; \
    else \
       echo "Delphi not found in downloads"; \
    fi

## Try installing XPLOR-NIH from downloads copy
RUN if  test -e downloads/xplor-nih-????-db.tar.gz \
      && test -e downloads/xplor-nih-????-Linux_x86_64.tar.gz; then \
        cd /opt ; \
        tar zxvf /app/downloads/xplor-nih-????-db.tar.gz; \
        tar zxvf /app/downloads/xplor-nih-????-Linux_x86_64.tar.gz; \
        cd xplor-nih*; \
        ./configure -symlinks /usr/local/bin; \
        rm /app/downloads/xplor-nih-????-db.tar.gz; \
        rm /app/downloads/xplor-nih-????-Linux_x86_64.tar.gz; \
        ## try to minimize xplorer installation
        rm -r python deprecated eginput tcl test; \ 
        cd /app ; \
        echo "XPlor NIH installed from downloads copy."; \
    else \
        echo "XPlor NIH not found in downloads."; \
    fi 

## duplicate, just in case requirements was updated without updating Dockerfile
RUN pip install -r requirements_extended.txt

ENV PYTHONPATH $PYTHONPATH:`pwd`

CMD ["python", "biskit/test.py"]

