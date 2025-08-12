FROM python:3.12-slim

#COPY requirements.txt .
#RUN pip3 install --no-cache-dir -r requirements.txt

RUN pip3 install --no-cache-dir \
    itk \
    itk-elastix \
    matplotlib \
    nibabel \
    numpy \
    pandas \
    scipy \
    scikit-image \
    scikit-learn
    
COPY . /app

ENTRYPOINT ["python3", "/app/fmriscript.py"]

CMD [""]