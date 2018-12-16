venv: setup.py
	python3 -m venv venv
	./venv/bin/pip install --upgrade pip
	./venv/bin/pip install -e .[tests]
	touch venv

test: venv
	./venv/bin/pytest --cov -rfsxEX --cov-report term-missing
