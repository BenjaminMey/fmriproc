FROM python:3.12-slim

COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

#RUN pip3 install --no-cache-dir \
#    itk \
#    itk-elastix \
#    matplotlib \
#    nibabel \
#    numpy \
#    pandas \
#    scipy \
#    scikit-image \
#    scikit-learn

ENV MPLCONFIGDIR=/tmp/matplotlib_config
RUN mkdir -p /tmp/matplotlib_config && chmod -R 777 /tmp/matplotlib_config

COPY . /app

WORKDIR /app

CMD ["python3", "/app/fmriscript.py", "-d", "/data", "-a", "SAG_T1_3D_MPRAGE.nii"]
