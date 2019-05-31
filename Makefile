all: test

install:
	pip install .

min-req-install:
	pip install -r min-requirements.txt

test:
	pip install pytest~=3.6.1 pytest-cov~=2.4
	pytest --cov=xena_gdc_etl
	pip freeze

lint:
	pip install flake8
	flake8
