FROM graik/biskitbase

RUN pip install --no-cache coverage coveralls

## Everything following ADD is not cached by Docker
## -> the following RUN commands have to be re-executed with every build
ADD . /app
WORKDIR /app

ENV AMBERHOME /opt/amber18
ENV PATH $PATH:$AMBERHOME/bin
ENV PYTHONPATH $PYTHONPATH:`pwd`

CMD ["python", "biskit/test.py -e fails old extra"]
