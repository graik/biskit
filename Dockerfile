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

CMD ["python", "biskit/test.py", "-e", "exe"]

