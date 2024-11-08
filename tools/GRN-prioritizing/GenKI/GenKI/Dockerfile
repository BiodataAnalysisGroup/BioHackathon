FROM python:3.9-slim-buster

WORKDIR GenKI/

COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt

COPY GenKI.py GenKI.py
CMD ["python3", "/GenKI.py"]
